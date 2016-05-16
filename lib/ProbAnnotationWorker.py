
import subprocess
import sys
import os
import shutil
import traceback
import time
import math
import re
import tempfile
from biop3.ProbModelSEED.ProbAnnotationParser import ProbAnnotationParser
from biokbase import log
from urllib2 import HTTPError
from ConfigParser import ConfigParser

# Exception thrown when no features are found in Genome object
class NoFeaturesError(Exception):
    pass

# Exception thrown when blast command failed
class BlastError(Exception):
    pass

# Exception thrown when there is an invalid number calculating likelihoods
class BadLikelihoodError(Exception):
    pass

# Exception thrown when a target id is not found in rolestring dictionary
class NoTargetIdError(Exception):
    pass

# Exception thrown when there are no gene IDs in Genome object
class NoGeneIdsError(Exception):
    pass

# Exception thrown when role not found in roleToTotalProb dictionary
class RoleNotFoundEror(Exception):
    pass

''' Worker that implements probabilistic annotation algorithm. '''

class ProbAnnotationWorker:

    def __init__(self, genomeId, context=None):
        ''' Initialize object.

            @param genomeId: Genome ID string for genome being annotated
            @param context: User context when used in a server
            @return Nothing
        '''

        # Save the genome ID (used for messages and temporary file names).
        self.genomeId = genomeId

        # Get the configuration variables.
        serviceName = os.environ.get('KB_SERVICE_NAME', 'ProbModelSEED')
        cfg = ConfigParser()
        cfg.read(os.path.join(os.environ.get('KB_TOP'), 'deployment.cfg'))
        self.config = dict()
        for nameval in cfg.items(serviceName):
            self.config[nameval[0]] = nameval[1]
        
        # Use the context from the server or build a context when used outside of a server.
        if context is not None:
            self.ctx = context
        else:
            self.ctx = dict()
            self.ctx['client_ip'] = '127.0.0.1'
            self.ctx['user_id'] = '-'
            self.ctx['module'] = serviceName
            self.ctx['method'] = '-'
            self.ctx['call_id'] = '-'

        # Create a logger.
        self.logger = log.log(serviceName, ip_address=True, authuser=True, module=True, method=True,
            call_id=True, logfile=self.config['mlog_log_file'])
        self.logger.set_log_level(int(self.config['mlog_log_level']))

        # Create a ProbAnnotationParser object for working with the static database files.
        self.dataParser = ProbAnnotationParser(self.config)

        # Get the static database files.  If the files do not exist and they are downloaded
        # from Shock, it can take a few minutes before they are ready.
        self.dataParser.getDatabaseFiles(self.logger, '')
        
        # Create a work directory for storing temporary files.
        if not os.path.exists(self.config['work_dir']):
            os.makedirs(self.config['work_dir'], 0775)
        self.workFolder = tempfile.mkdtemp(dir=self.config['work_dir'], prefix='')

        return

    def genomeToFasta(self, features):

        ''' Convert the features from a genome into an amino-acid FASTA file (for BLAST purposes).

            @param features: List of features with protein sequences
            @return Path to fasta file with query proteins
            @raise NoFeaturesError when list of features is empty
        '''

        # Make sure the genome has features.
        if len(features) == 0:
            raise NoFeaturesError('Genome %s has no features. Did you forget to run annotate_genome?\n' %(self.genomeId))
    
        # Run the list of features to build the fasta file.
        self._log(log.DEBUG, 'Creating protein fasta file for genome '+self.genomeId)
        fastaFile = os.path.join(self.workFolder, '%s.faa' %(self.genomeId))
        with open(fastaFile, 'w') as handle:
            numProteins = 0
            for feature in features:
                # Not a protein-encoding gene
                if 'protein_translation' not in feature:
                    continue
                handle.write('>%s\n%s\n' %(feature['id'], feature['protein_translation']))
                numProteins += 1
        
        self._log(log.DEBUG, 'Wrote %d protein sequences to "%s"' %(numProteins, fastaFile))
        return fastaFile
        
    def runBlast(self, queryFile):

        ''' A simplistic wrapper to search for the query proteins against the subsystem proteins.

            @param queryFile: Path to fasta file with query proteins
            @return Path to output file from search program
            @raise BlastError when there is a problem running the search program
        '''

        # Generate path to output file.  Output format 6 is tab-delimited format.
        blastResultFile = os.path.join(self.workFolder, '%s.blastout' %(self.genomeId))

        # Build the command based on the configured search program.
        if self.config['search_program'] == 'usearch':
            args = [ self.config['search_program_path'], '-ublast', queryFile,
                     '-db', self.dataParser.SearchFiles['protein_udb_file'],
                     '-evalue', self.config['search_program_evalue'],
                     '-accel', self.config['usearch_accel'],
                     '-threads', self.config['search_program_threads'],
                     '-blast6out', blastResultFile ]
        else:
            args = [ self.config['search_program_path'], '-query', queryFile,
                     '-db', self.dataParser.DataFiles['protein_fasta_file'],
                     '-outfmt', '6', '-evalue', self.config['search_program_evalue'],
                     '-num_threads', self.config['search_program_threads'],
                     '-out', blastResultFile ]

        # Run the command to search for proteins against subsystem proteins.
        cmd = ' '.join(args)
        self._log(log.DEBUG, 'Started protein search with command: '+cmd)
        try:
            proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            (stdout, stderr) = proc.communicate()
            if proc.returncode < 0:
                message = '"%s" was terminated by signal %d' %(args[0], -proc.returncode)
                raise BlastError(message)
            else:
                if proc.returncode > 0:
                    details = '"%s" failed with return code %d\nCommand: "%s"\nStdout: "%s"\nStderr: "%s"' \
                        %(args[0], proc.returncode, cmd, stdout, stderr)
                    raise BlastError(details)
        except OSError as e:
            message = 'Failed to run "%s": %s' %(args[0], e.strerror)
            raise BlastError(message)
        self._log(log.DEBUG, 'Finished protein search')

        return blastResultFile
    
    def rolesetProbabilitiesMarble(self, blastResultFile):

        ''' Calculate the probabilities of rolesets from the BLAST results.

            A roleset is each possible combination of roles implied by the functions
            of the proteins in subsystems.  The output is a dictionary keyed by
            query gene of lists of tuples where each tuple contains (1) roleset
            string, and (2) likelihood value.  The roleset string is a concatenation
            of all of the roles of a protein with a single function (order does
            not matter).
    
            @param blastResultFile: Path to output file from BLAST
            @return Dictionary keyed by query gene of list of tuples with roleset and likelihood
            @raise BadLikelihoodError when there is math error calculating a likelihood
            @raise NoTargetIdError when target ID in search results is not found in rolestrings
        '''

        self._log(log.DEBUG, 'Started marble-picking on rolesets for genome '+self.genomeId)
    
        # Read in the target roles (this function returns the roles as lists!)
        targetIdToRole, targetRoleToId = self.dataParser.readFidRoleFile(self.dataParser.DataFiles['otu_fid_role_file'])
    
        # Convert the lists of roles into "rolestrings" (sort the list so that order doesn't matter)
        # in order to deal with the case where some of the hits are multi-functional and others only have
        # a single function.
        targetIdToRoleString = dict()
        for target in targetIdToRole:
            stri = self.config['separator'].join(sorted(targetIdToRole[target]))
            targetIdToRoleString[target] = stri

        # Parse the output from BLAST which returns a dictionary keyed by query gene of a list
        # of tuples with target gene and score.
        # query --> [ (target1, score 1), (target 2, score 2), ... ]
        idToTargetList = self.dataParser.parseBlastOutput(blastResultFile)
    
        # This is a holder for all of our results which is a dictionary keyed by query gene
        # of a list of tuples with roleset and likelihood.
        # query -> [ (roleset1, likelihood_1), (roleset2, likelihood_2), ...]
        rolestringTuples = dict()

        # For each query gene we calculate the likelihood of each possible rolestring
        # See equation 2 in the paper ("Calculating annotation likelihoods" section).
        for query in idToTargetList:
            # First we need to know the maximum score for this gene.
            # I have no idea why but I'm pretty sure Python is silently turning the second
            # element of these tuples into strings.  That's why I turn them back to floats.
            maxscore = 0
            for tup in idToTargetList[query]:
                if float(tup[1]) > maxscore:
                    maxscore = float(tup[1])
    
            # Now we calculate the cumulative squared scores for each possible rolestring.
            # This along with pseudocount*maxscore is equivalent to multiplying all scores
            # by themselves and then dividing by the max score.
            # This is done to avoid some pathological cases and give more weight to higher-scoring hits
            # and not let much lower-scoring hits \ noise drown them out.
            # Build a dictionary keyed by rolestring of the sum of squares of the log-scores.
            rolestringToScore = dict()
            for tup in idToTargetList[query]:
                try:
                    rolestring = targetIdToRoleString[tup[0]]
                except KeyError:
                    message = 'Target id %s from search results file had no roles in rolestring dictionary' %(tup[0])
                    raise NoTargetIdError(message)
                if rolestring in rolestringToScore:
                    rolestringToScore[rolestring] += (float(tup[1]) ** 2)
                else:
                    rolestringToScore[rolestring] = (float(tup[1]) ** 2)
    
            # Calculate the likelihood that this gene has the given functional annotation.
            # Start with the denominator which is the sum of squares of the log-scores for
            # all possible rolestrings.
            denom = float(self.config['pseudo_count']) * maxscore
            for stri in rolestringToScore:
                denom += rolestringToScore[stri]
            if math.isnan(denom):
                message = 'Denominator in likelihood calculation for gene %s is NaN %f' %(query, denom)
                raise BadLikelihoodError(message)

            # The numerators are the sum of squares for each rolestring.
            # Calculate the likelihood for each rolestring and store in the output dictionary.
            for stri in rolestringToScore:
                p = rolestringToScore[stri] / denom
                if math.isnan(p):
                    message = 'Likelihood for rolestring %s in gene %s is NaN based on score %f' %(stri, query, rolestringToScore[stri])
                    raise BadLikelihoodError(message)
                if query in rolestringTuples:
                    rolestringTuples[query].append( (stri, p) )
                else:
                    rolestringTuples[query] = [ (stri, p) ]
    
        # Save the generated data when debug is turned on.
        if self.logger.get_log_level() >= log.DEBUG2:
            rolesetProbabilityFile = os.path.join(self.workFolder, '%s.rolesetprobs' %(self.genomeId))
            with open(rolesetProbabilityFile, 'w') as handle:
                for query in rolestringTuples:
                    for tup in rolestringTuples[query]:
                        handle.write('%s\t%s\t%1.4f\n' %(query, tup[0], tup[1]))
            
        self._log(log.DEBUG, 'Finished marble-picking on %d rolesets for genome %s' %(len(rolestringTuples), self.genomeId))
        return rolestringTuples
            
    def rolesetProbabilitiesToRoleProbabilities(self, queryToTuplist):
        ''' Compute probability of each role from the rolesets for each query protein.

            At the moment the strategy is to take any set of rolestrings containing
            the same roles and add their probabilities.  So if we have hits to both
            a bifunctional enzyme with R1 and R2, and hits to a monofunctional enzyme
            with only R1, R1 ends up with a greater probability than R2.

            I had tried to normalize to the previous sum but I need to be more careful
            than that (I'll put it on my TODO list) because if you have e.g. one hit
            to R1R2 and one hit to R3 then the probability of R1 and R2 will be unfairly
            brought down due to the normalization scheme.

            @param queryToTuplist: Dictionary keyed by query gene of list of tuples with roleset and likelihood
            @return List of tuples with query gene, role, and likelihood
        '''

        self._log(log.DEBUG, 'Started computing role probabilities for genome '+self.genomeId)

        # Start with an empty list.
        roleProbs = list()

        # Iterate over all of the query genes in the dictionary.
        # querygene -> [ (roleset1, likelihood_1), (roleset2, likelihood_2), ...]
        for query in queryToTuplist:
            # This section actually does the conversion of likelihoods.
            # See equation 3 in the paper ("Calculating reaction likelihoods" section).
            queryRolesToProbs = dict()
            for tup in queryToTuplist[query]:
                rolelist = tup[0].split(self.config['separator'])
                # Add up all the instances of each particular role on the list.
                for role in rolelist:
                    if role in queryRolesToProbs:
                        queryRolesToProbs[role] += tup[1]
                    else:
                        queryRolesToProbs[role] = tup[1]

            # Add them to the array.
            for role in queryRolesToProbs:
                roleProbs.append( (query, role, queryRolesToProbs[role]) )

        # Save the generated data when debug is turned on.
        if self.logger.get_log_level() >= log.DEBUG2:
            role_probability_file = os.path.join(self.workFolder, '%s.roleprobs' %(self.genomeId))
            with open(role_probability_file, "w") as handle:
                for tuple in roleProbs:
                    handle.write('%s\t%s\t%s\n' %(tuple[0], tuple[1], tuple[2]))

        self._log(log.DEBUG, 'Finished computing %d role probabilities for genome %s' %(len(roleProbs), self.genomeId))

        return roleProbs

    def totalRoleProbabilities(self, roleProbs):
        ''' Given the likelihood that each gene has each role, estimate the likelihood
            that the entire ORGANISM has that role.

            To avoid exploding the likelihoods with noise, I just take the maximum
            likelihood of any query gene having a function and use that as the
            likelihood that the function exists in the cell.

            A gene is assigned to a role if it is within DILUTION_PERCENT of the maximum
            probability. DILUTION_PERCENT can be adjusted in the config file. For each
            role the maximum likelihood and the estimated set of genes that perform that
            role are linked with an OR relationship to form a Boolean Gene-Function
            relationship.

            @param roleProbs List of tuples with query gene, role, and likelihood
            @return List of tuples with role, likelihood, and estimated set of genes that perform the role
            @raise RoleNotFoundError when role is not placed properly in roleToTotalProb dictionary
        '''

        self._log(log.DEBUG, 'Started generating whole-cell role probabilities for genome '+self.genomeId)

        # Find maximum likelihood among all query genes for each role.
        # This is assumed to be the likelihood of that role occurring in the organism as a whole.
        roleToTotalProb = dict()
        for tuple in roleProbs:
            if tuple[1] in roleToTotalProb:
                if float(tuple[2]) > roleToTotalProb[tuple[1]]:
                    roleToTotalProb[tuple[1]] = float(tuple[2])
            else:
                roleToTotalProb[tuple[1]] = float(tuple[2])

        # Get the genes within DILUTION_PERCENT percent of the maximum
        # likelihood and assert that these are the most likely genes responsible for that role.
        # (note - DILUTION_PERCENT is defined in the config file)
        # This produces a dictionary from role to a list of genes
        # See equation 4 in the paper ("Calculating reaction likelihoods" section).
        roleToGeneList = dict()
        for tuple in roleProbs:
            if tuple[1] not in roleToTotalProb:
                message = 'Role %s not placed properly in roleToTotalProb dictionary?' %(tuple[1])
                self._log(log.ERR, message)
                raise RoleNotFoundError(message)
            if float(tuple[2]) >= float(self.config['dilution_percent'])/100.0 * roleToTotalProb[tuple[1]]:
                if tuple[1] in roleToGeneList:
                    roleToGeneList[tuple[1]].append(tuple[0])
                else:
                    roleToGeneList[tuple[1]] = [ tuple[0] ]

        # Build the array of total role probabilities.
        totalRoleProbs = list()
        for role in roleToTotalProb:
            gpr = ' or '.join(list(set(roleToGeneList[role])))
            # We only need to group these if there is more than one of them (avoids extra parenthesis when computing complexes)
            if len(list(set(roleToGeneList[role]))) > 1:
                gpr = '(' + gpr + ')'
            totalRoleProbs.append( (role, roleToTotalProb[role], gpr ) )

        # Save the generated data when debug is turned on.
        if self.logger.get_log_level() >= log.DEBUG2:
            total_role_probability_file = os.path.join(self.workFolder, '%s.cellroleprob' %(self.genomeId))
            with open(total_role_probability_file, "w") as handle:
                for tuple in totalRoleProbs:
                    handle.write('%s\t%s\t%s\n' %(tuple[0], tuple[1], tuple[2]))

        self._log(log.DEBUG, 'Finished generating %d whole-cell role probabilities for genome %s' %(len(totalRoleProbs), self.genomeId))

        return totalRoleProbs

    def complexProbabilities(self, totalRoleProbs, complexesToRequiredRoles = None):
        ''' Compute the likelihood of each protein complex from the likelihood of each role.

            A protein complex represents a set functional roles that must all be present
            for a complex to exist.  The likelihood of the existence of a complex is
            computed as the minimum likelihood of the roles within that complex (ignoring
            roles not represented in the subsystems).

            For each protein complex, the likelihood, type, list of roles not in the
            organism, and list of roles not in subsystems is returned.  The type is a
            string with one of the following values:

            CPLX_FULL - All roles found in organism and utilized in the complex
            CPLX_PARTIAL - Only some roles found in organism and only those roles that
                were found were utilized. Note this does not distinguish between not
                there and not represented for roles that were not found
            CPLX_NOTTHERE - Likelihood is 0 because the genes aren't there for any of
                the subunits
            CPLX_NOREPS - Likelihood is 0 because there are no representative genes in
                the subsystems for any of the subunits
            CPLX_NOREPS_AND_NOTTHERE - Likelihood is 0 because some genes aren't there
                for any of the subunits and some genes have no representatives

            @param totalRoleProbs: List of tuples with role, likelihood, and estimated set
                of genes that perform the role
            @param complexesToRequiredRoles: Dictionary keyed by complex ID to the roles
                involved in forming that complex. If it is None we read it from the CDMI
                files we downloaded, otherwise we use the provided dictionary.
            @return List of tuples with complex ID, likelihood, type, list of roles not in
                organism, list of roles not in subsystems, and boolean Gene-Protein
                relationship
        '''

        self._log(log.DEBUG, 'Started computing complex probabilities for '+self.genomeId)

        # Get the mapping from complexes to roles if it isn't already provided.
        if complexesToRequiredRoles is None:
            complexesToRequiredRoles, rolesToComplexes = self.dataParser.readComplexRoleFile(self.dataParser.DataFiles['complex_role_file'])
            self._log(log.DEBUG, 'Found %d complex to role mappings in %s' %(len(complexesToRequiredRoles), self.dataParser.DataFiles['complex_role_file']))

        # Get the subsystem roles (used to distinguish between NOTTHERE and NOREPS).
        otu_fidsToRoles, otu_rolesToFids = self.dataParser.readFidRoleFile(self.dataParser.DataFiles['otu_fid_role_file'])
        allroles = set()
        for fid in otu_fidsToRoles:
            for role in otu_fidsToRoles[fid]:
                allroles.add(role)

        # Build two dictionaries, both keyed by role, one mapping the role to its
        # likelihood and one mapping to the gene list.
        rolesToProbabilities = dict()
        rolesToGeneList = dict()
        for tuple in totalRoleProbs:
            rolesToProbabilities[tuple[0]] = float(tuple[1]) # can skip the float()?
            rolesToGeneList[tuple[0]] = tuple[2]

        # Iterate over complexes and compute complex probabilities from role probabilities.
        # Separate out cases where no genes seem to exist in the organism for the reaction
        # from cases where there is a database deficiency.
        # See equation 5 in the paper ("Calculating reaction likelihoods" section).
        SEPARATOR = self.config['separator']
        complexProbs = list()
        for cplx in complexesToRequiredRoles:
            allCplxRoles = complexesToRequiredRoles[cplx]
            availRoles = list() # Roles that may have representatives in the query organism
            unavailRoles = list() # Roles that have representatives but that are not apparently in the query organism
            noexistRoles = list() # Roles with no representatives in the subsystems
            for role in complexesToRequiredRoles[cplx]:
                if role not in allroles:
                    noexistRoles.append(role)
                elif role not in rolesToProbabilities:
                    unavailRoles.append(role)
                else:
                    availRoles.append(role)
            TYPE = ""
            GPR = ""
            if len(noexistRoles) == len(allCplxRoles):
                TYPE = "CPLX_NOREPS"
                complexProbs.append( (cplx, 0.0, TYPE, self.config["separator"].join(unavailRoles), self.config["separator"].join(noexistRoles), GPR) )
                continue
            if len(unavailRoles) == len(allCplxRoles):
                TYPE = "CPLX_NOTTHERE"
                complexProbs.append( (cplx, 0.0, TYPE, self.config["separator"].join(unavailRoles), self.config["separator"].join(noexistRoles), GPR) )
                continue
            # Some had no representatives and the rest were not found in the cell
            if len(unavailRoles) + len(noexistRoles) == len(allCplxRoles):
                TYPE = "CPLX_NOREPS_AND_NOTTHERE"
                complexProbs.append( (cplx, 0.0, TYPE, self.config["separator"].join(unavailRoles), self.config["separator"].join(noexistRoles), GPR) )
                continue
            # Otherwise at least one of them is available
            if len(availRoles) == len(allCplxRoles):
                TYPE = "CPLX_FULL"
            elif len(availRoles) < len(allCplxRoles):
                TYPE = "CPLX_PARTIAL_%d_of_%d" %(len(availRoles), len(allCplxRoles))

            # Link individual functions in complex with an AND relationship to form a
            # Boolean Gene-Protein relationship.
#            partialGprList = [ "(" + s + ")" for s in [ rolesToGeneList[f] for f in availRoles ] ]
            partialGprList = [ rolesToGeneList[f] for f in availRoles ]
            GPR = " and ".join( list(set(partialGprList)) )

            if GPR != "" and len(list(set(partialGprList))) > 1:
                GPR = "(" + GPR + ")"

            # Find the minimum probability of the different available roles (ignoring ones
            # that are apparently missing) and call that the complex likelihood.
            minp = 1000
            for role in availRoles:
                if rolesToProbabilities[role] < minp:
                    minp = rolesToProbabilities[role]
            complexProbs.append( (cplx, minp, TYPE, self.config["separator"].join(unavailRoles), self.config["separator"].join(noexistRoles), GPR) )

        # Save the generated data when debug is turned on.
        if self.logger.get_log_level() >= log.DEBUG2:
            complex_probability_file = os.path.join(self.workFolder, "%s.complexprob" %(self.genomeId))
            with open(complex_probability_file, "w") as handle:
                for tuple in complexProbs:
                    handle.write("%s\t%1.4f\t%s\t%s\t%s\t%s\n" %(tuple[0], tuple[1], tuple[2], tuple[3], tuple[4], tuple[5]))

        self._log(log.DEBUG, 'Finished computing complex probabilities for '+self.genomeId)
        return complexProbs

    def reactionProbabilities(self, complexProbs, rxnsToComplexes = None):
        ''' Estimate the likelihood of reactions from the likelihood of complexes.

            The reaction likelihood is computed as the maximum likelihood of complexes
            that perform that reaction.

            If the reaction has no complexes it won't even be in this file because of the way
            I set up the call... I could probably change this so that I get a list of ALL reactions
            and make it easier to catch issues with reaction --> complex links in the database.
            Some of the infrastructure is already there (with the TYPE).

            @param complexProbs: List of tuples with complex ID, likelihood, type, list of
                roles not in organism, list of roles not in subsystems, and boolean
                Gene-Protein relationship
            @param rxnsToComplexes: Dictionary keyed by reaction ID to a list of catalyzing
                complexes. If it is None we read it from the CDMI files we downloaded,
                otherwise we use the provided dictionary.
            @return List of tuples with reaction ID, likelihood, reaction type, complex info,
                and gene-protein-reaction relationship
        '''

        self._log(log.DEBUG, 'Started computing reaction probabilities for '+self.genomeId)

        # Build a dictionary keyed by complex ID of tuples with likelihood, type, and GPR.
        # Note we don't need to use the list of roles not in organism and list of roles
        # not in subsystems.
        # cplx --> {likelihood, type, GPR}
        cplxToTuple = dict()
        for tuple in complexProbs:
            cplxToTuple[tuple[0]] = ( tuple[1], tuple[2], tuple[5] )
        
        # Get the mapping from reactions to complexes if it isn't already provided.
        if rxnsToComplexes is None:
            rxnsToComplexes = self.dataParser.readReactionComplexFile(self.dataParser.DataFiles['reaction_complex_file'])

        # Take the MAXIMUM likelihood of complexes catalyzing a particular reaction
        # and call that the reaction likelihood.
        # See equation 6 in the paper ("Calculating reaction likelihoods" section).
        reactionProbs = list()
        for rxn in rxnsToComplexes:
            TYPE = "NOCOMPLEXES"
            rxnComplexes = rxnsToComplexes[rxn]
            maxProb = 0
            GPR = ""
            complexList = list()
            for cplx in rxnComplexes:
                if cplx in cplxToTuple:
                    # Complex1 (P1; TYPE1) ///Complex2 (P2; TYPE2) ...
                    complexList.append( [ cplx, cplxToTuple[cplx][0], cplxToTuple[cplx][1] ])
                    TYPE = 'HASCOMPLEXES'
            complexString = ''
            if len(complexList) > 0:
                complexList.sort(key=lambda tup: tup[1], reverse=True)
                maxProb = complexList[0][1]
                for complex in complexList:
                    complexString += '%s (%1.4f; %s)%s' %(complex[0], complex[1], complex[2], self.config['separator'])
                complexString = complexString[:-len(self.config['separator'])] # Remove the final separator

            # Iterate separately to get a GPR. We want to apply a cutoff here too to avoid
            # a complex with 80% probability being linked by OR to another with a 5%
            # probability.  For now I've implemented using the same cutoff as we used for
            # which genes go with a role.
            cplxGprs = []
            for cplx in rxnComplexes:
                if cplx in cplxToTuple:
                    if cplxToTuple[cplx][0] < maxProb * float(self.config["dilution_percent"])/100.0:
                        continue
                    cplxGprs.append(cplxToTuple[cplx][2])
            if len(cplxGprs) > 0:
                GPR = " or ".join( list(set(cplxGprs)) )

            # Add everything to the final list.
            reactionProbs.append( [rxn, maxProb, TYPE, complexString, GPR] )

        # Save the generated data when debug is turned on.
        if self.logger.get_log_level() >= log.DEBUG2:
            reaction_probability_file = os.path.join(self.workFolder, "%s.rxnprobs" %(self.genomeId))
            with open(reaction_probability_file, "w") as handle:
                for tuple in reactionProbs:
                    handle.write("%s\t%1.6f\t%s\t%s\t%s\n" %(tuple[0], tuple[1], tuple[2], tuple[3], tuple[4]))

        self._log(log.DEBUG, 'Finished computing reaction probabilities for '+self.genomeId)
        return reactionProbs

    def cleanup(self):
        ''' Cleanup the work folder.

            @return Nothing
        '''

        shutil.rmtree(self.workFolder)
        return

    def _log(self, level, message):
        ''' Log a message to the system log.

            @param level: Message level (INFO, WARNING, etc.)
            @param message: Message text
            @return Nothing
        '''

        # Log the message.
        self.logger.log_message(level, message, self.ctx['client_ip'], self.ctx['user_id'], self.ctx['module'],
                                self.ctx['method'], self.ctx['call_id'])
        return
