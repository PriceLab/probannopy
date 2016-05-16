
# Read and write data files
import os
import sys
import math
import subprocess
import json
import traceback
import time
from shock import Client as ShockClient
from biokbase import log

# E values of less than 1E-200 are treated as 1E-200 to avoid log of 0 issues.
MIN_EVALUE = 1E-200

# Exception thrown when there is an invalid sources configuration variable
class BadSourceError(Exception):
    pass

# Exception thrown when makeblastdb command failed
class MakeblastdbError(Exception):
    pass

# Exception thrown when static database file is missing from Shock.
class MissingFileError(Exception):
    pass

# Exception thrown when static database files are not ready
class NotReadyError(Exception):
    pass

''' Read and write data files. '''

class ProbAnnotationParser:
    
    def __init__(self, config):
        ''' Initialize the object.
        
            @param config Dictionary of configuration variables
        '''

        # Save the configuration variables related to data files.
        self.dataFolderPath = config['data_dir']
        self.separator = config['separator']
        self.searchProgram = config['search_program']
        self.searchProgramPath = config['search_program_path']
        self.shockURL = config['shock_url']
        self.loadDataOption = config['load_data_option']

        # Create a dictionary with the valid sources and initialize to not set.
        self.sources = { 'cdm': dict(), 'kegg': dict() }

        # The list of sources are separated with a colon.
        sourceList = config['data_sources'].split(':')
        if len(sourceList) == 0:
            raise BadSourceError('List of data sources is empty')
        
        # Check for valid sources and mark each one that is in the list.
        for src in sourceList:
            if src not in self.sources:
                raise BadSourceError('Source "%s" is not supported' %(src))
            if src == 'cdm':
                self.sources['cdm']['subsystem_fid_file'] = os.path.join(self.dataFolderPath, 'CDM_SUBSYSTEM_FID')
                self.sources['cdm']['dlit_fid_file'] = os.path.join(self.dataFolderPath, 'CDM_DLIT_FID')
                self.sources['cdm']['concatenated_fid_file'] = os.path.join(self.dataFolderPath, 'CDM_ALL_FID')
                self.sources['cdm']['fid_role_file'] = os.path.join(self.dataFolderPath, 'CDM_ALL_FID_ROLE')
                self.sources['cdm']['otu_id_file'] = os.path.join(self.dataFolderPath, 'OTU_GENOME_IDS')
                self.sources['cdm']['otu_fid_role_file'] = os.path.join(self.dataFolderPath, 'CDM_OTU_FID_ROLE')
                self.sources['cdm']['protein_fasta_file'] = os.path.join(self.dataFolderPath, 'CDM_PROTEIN_FASTA')
                self.sources['cdm']['complex_role_file'] = os.path.join(self.dataFolderPath, 'CDM_COMPLEX_ROLE')
                self.sources['cdm']['reaction_complex_file'] = os.path.join(self.dataFolderPath, 'CDM_REACTION_COMPLEX')
            if src == 'kegg':
                self.sources['kegg']['reaction_file'] = os.path.join(self.dataFolderPath, 'KEGG_KBASE_REACTION')
                self.sources['kegg']['otu_fid_role_file'] = os.path.join(self.dataFolderPath, 'KEGG_OTU_FID_ROLE')
                self.sources['kegg']['protein_fasta_file'] = os.path.join(self.dataFolderPath, 'KEGG_PROTEIN_FASTA')
                self.sources['kegg']['complex_role_file'] = os.path.join(self.dataFolderPath, 'KEGG_COMPLEX_ROLE')
                self.sources['kegg']['reaction_complex_file'] = os.path.join(self.dataFolderPath, 'KEGG_REACTION_COMPLEX')

        # One of the source must be central data model.
        if len(self.sources['cdm']) == 0:
            raise BadSourceError('One of the data sources must be cdm')

        # Paths to files for tracking status of static database files.
        self.StatusFiles = dict()
        self.StatusFiles['status_file'] = os.path.join(self.dataFolderPath, 'staticdata.status')
        self.StatusFiles['cache_file'] = os.path.join(self.dataFolderPath, 'staticdata.cache')

        # Paths to files with source data.
        self.DataFiles = dict()
        self.DataFiles['otu_fid_role_file'] = os.path.join(self.dataFolderPath, 'OTU_FID_ROLE')
        self.DataFiles['protein_fasta_file'] = os.path.join(self.dataFolderPath, 'PROTEIN_FASTA')
#         self.DataFiles['complex_role_file'] = os.path.join(self.dataFolderPath, 'COMPLEX_ROLE')
#         self.DataFiles['reaction_complex_file'] = os.path.join(self.dataFolderPath, 'REACTION_COMPLEX')

        # Paths to files for searching for proteins.
        self.SearchFiles = dict()
        if self.searchProgram == 'usearch':
            self.SearchFiles['protein_udb_file'] = os.path.join(self.dataFolderPath, 'PROTEIN.udb')
        else:
            self.SearchFiles['protein_otu_index_file'] = os.path.join(self.dataFolderPath, 'PROTEIN_FASTA.pin')
            self.SearchFiles['protein_otu_sequence_file'] = os.path.join(self.dataFolderPath, 'PROTEIN_FASTA.psq')
            self.SearchFiles['protein_otu_header_file'] = os.path.join(self.dataFolderPath, 'PROTEIN_FASTA.phr')

        # Create the data folder if it does not exist.
        if not os.path.exists(config['data_dir']):
            os.makedirs(config['data_dir'], 0775)

        return

    # The OTU ID file is a list of representative OTU genome IDs.  Each line has these fields:
    #   1. Genome ID in KBase format (e.g. kb|g.0)
    #   2. Flag indicating if the genome is a prokaryote (1 means yes, 0 means no)
    
    def readOtuData(self, filename):
        ''' Read data from the representative OTU genome ID file.

            @param filename: Path to OTU genome ID file
            @return List of all OTU genome IDs, list of prokaryote OTU genome IDs
        '''
    
        otus = list()
        prokotus = list()
        with open(filename, 'r') as handle:
            for line in handle:
                fields = line.strip('\r\n').split('\t')
                otus.append(fields[0])
                if int(fields[1]) == 1:
                    prokotus.append(fields[0])
        return otus, prokotus
    
    def writeOtuData(self, filename, otus, prokotus):
        ''' Write data to the representative OTU genome ID file.

            @param filename: Path to OTU genome ID file
            @param otus: List of all OTU genome IDs
            @param prokotus: List of prokaryote OTU genome IDs
            @return Nothing
        '''

        with open(filename, 'w') as handle:
            for otu in otus:
                if otu in prokotus:
                    handle.write('%s\t%d\n' %(otu, 1))
                else:
                    handle.write('%s\t%d\n' %(otu, 0))
        return
    
    # A feature ID file is a list of features IDs.  Each line has one field that is
    # the feature ID in KBase format (e.g. kb|g.3.peg.541).

    def readFeatureIdFile(self, filename):
        ''' Read data from a feature ID file.
        
            @param filename: Path to feature ID file
            @return List of feature IDs.
        '''
        featureIds = list()
        with open(filename, 'r') as handle:
            for line in handle:
                featureIds.append(line.strip('\r\n'))
        return featureIds

    def writeFeatureIdFile(self, filename, featureIds):
        ''' Write data to a feature ID file.
        
            @param filename: Path to feature ID file
            @param featureIds: List of feature IDs
            @return Nothing
        '''
        
        with open(filename, 'w') as handle:
            for index in range(len(featureIds)):
                handle.write('%s\n' %(featureIds[index]))
        return

    # A feature ID to role file is a mapping of feature IDs to functional roles.
    # Each line has these fields:
    #   1. Feature ID in KBase format (e.g. kb|g.0.peg.2094)
    #   2. List of names of functional roles (e.g. Conserved ATP-binding protein YghS)
    #
    # Note that functional roles must be separated by a string that does not occur in any role.
    
    def readFidRoleFile(self, filename):
        ''' Read data from a feature ID to role file.
        
            @param filename: Path to feature ID to role file
            @return Dictionary mapping a feature ID to list of names of roles,
                dictionary mapping a role to feature ID
        '''
        
        fidsToRoles = dict()
        with open(filename, 'r') as handle:
            for line in handle:
                fields = line.strip('\r\n').split('\t')
                roles = fields[1].split(self.separator)
                if fields[0] in fidsToRoles:
                    fidsToRoles[fields[0]] += roles
                else:
                    fidsToRoles[fields[0]] = roles

        rolesToFids = dict()
        for fid in fidsToRoles:
            roles = fidsToRoles[fid]
            for role in roles:
                if role in rolesToFids:
                    rolesToFids[role].append(fid)
                else:
                    rolesToFids[role] = [ fid ]
    
        return fidsToRoles, rolesToFids

    def writeFidRoleFile(self, filename, fidsToRoles):
        ''' Write data to a feature ID to role file.
        
            @param filename: Path to feature ID to role file
            @param fidsToRoles: Dictionary mapping a feature ID to list of names of roles
            @return Nothing
        '''

        with open(filename, 'w') as handle:
            for fid in fidsToRoles:
                try:
                    handle.write('%s\t%s\n' %(fid, self.separator.join(fidsToRoles[fid])))
                except UnicodeEncodeError:
                    roles = list()
                    for index in range(len(fidsToRoles[fid])):
                        roles.append(fidsToRoles[fid][index].encode('ascii', 'replace'))
                    handle.write('%s\t%s\n' %(fid, self.separator.join(roles)))
        return

    # A protein FASTA file contains the amino acid sequences for a set of feature IDs.
    
    def writeProteinFastaFile(self, filename, fidsToSeqs):
        ''' Write data to a protein FASTA file.
        
            @param filename: Path to protein FASTA file
            @param fidsToSeqs: Dictionary mapping a feature ID to amino acid sequence
            @return Nothing
        '''
    
        # Sort the fids so that fasta files containing the same proteins hash to the same MD5 (for
        # data provenance purposes)
        with open(filename, 'w') as handle:
            for fid in sorted(fidsToSeqs.keys()):
                handle.write('>%s\n%s\n' %(fid, fidsToSeqs[fid]))
        return
    
    def buildSearchDatabase(self):
        ''' Build a search database for the configured search program.

            @note Make sure the subsystem FASTA file is available.
            @return Nothing
        '''

        # Build the command based on the configured search program.
        if self.searchProgram == 'usearch':
            args = [ self.searchProgramPath, '-makeudb_ublast', self.DataFiles['protein_fasta_file'], '-output', self.SearchFiles['protein_udb_file'] ]
        else:
            args = [ '/usr/bin/makeblastdb', '-in', self.DataFiles['protein_fasta_file'], '-dbtype', 'prot' ]
    
        # Run the command to compile the database from the subsystem fasta file.
        try:
            proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            (stdout, stderr) = proc.communicate()
            if proc.returncode < 0:
                cmd = ' '.join(args)
                raise MakeblastdbError('"%s" was terminated by signal %d' %(cmd, -proc.returncode))
            else:
                if proc.returncode > 0:
                    cmd = ' '.join(args)
                    details = '"%s" failed with return code %d:\nCommand: "%s"\nStdout: "%s"\nStderr: "%s"' \
                        %(args[0], proc.returncode, cmd, stdout, stderr)
                    raise MakeblastdbError(details)
        except OSError as e:
            cmd = ' '.join(args)
            raise MakeblastdbError('Failed to run "%s": %s' %(cmd, e.strerror))
        return
    
    def parseBlastOutput(self, blastResultsPath):
        ''' Read BLAST results file and store in a convenient structure.

            The results file is in BLAST output format 6 where each line describes an alignment
            found by the search program.  A line has 12 tab delimited fields: (1) query label,
            (2) target label, (3) percent identity, (4) alignment length, (5) number of
            mismatches, (6) number of gap opens, (7) 1-based position of start in query,
            (8) 1-based position of end in query, (9) 1-based position of start in target,
            (10) 1-based position of end in target, (11) e-value, and (12) bit score.
 
            @note Score is the negative log E-value
            @param blastResultsPath Path to BLAST results file
            @return Dictionary mapping query ID to tuple of target ID and score
        '''
    
        idToTargetList = dict()
        for line in open(blastResultsPath, 'r'):
            fields = line.strip('\r\n').split('\t')
            queryid = fields[0]
            targetid = fields[1]
            if float(fields[11]) < 0.0: # Throw out alignments with a negative bit score
                print 'throwing out %s' %(line)
                continue
            logeval = -1.0 * math.log10(float(fields[10]) + MIN_EVALUE)
            tup = ( targetid, logeval )
            if queryid in idToTargetList:
                idToTargetList[queryid].append( tup )
            else:
                idToTargetList[queryid] = [ tup ]
        return idToTargetList
    
    # A complexes to roles file contains a mapping of complex IDs to functional roles.
    # Each line has these fields:
    #   1. Complex ID in KBase format (e.g. kb|cpx.1048)
    #   2. List of functional roles for the complex
    #
    # Note that functional role names must be separated by a string that does not occur in any name.

    def readComplexRoleFile(self, filename):
        ''' Read data from a complex to role file.

            @param filename: Path to complex to roles file
            @return Dictionary mapping a complex ID to list of names of functional roles,
                dictionary mapping a role to a list of complex IDs
        '''

        complexToRoles = dict()
        roleToComplexes = dict()
        with open(filename, 'r') as handle:
            for line in handle:
                fields = line.strip('\r\n').split('\t')
                complex = fields[0]
                roles = fields[1].split(self.separator)
                if complex not in complexToRoles:
                    complexToRoles[complex] = roles
                else:
                    complexToRoles[complex] += roles
                for rindex in range(len(roles)):
                    if roles[rindex] not in roleToComplexes:
                        roleToComplexes[roles[rindex]] = list()
                    roleToComplexes[roles[rindex]].append(complex)
        return complexToRoles, roleToComplexes

    def writeComplexRoleFile(self, filename, complexToRoles):
        ''' Write data to a complex to role file.

            @param filename: Path to complex to role file
            @param complexToRoles: Dictionary mapping a complex ID to list of names of functional roles
            @return Nothing
        '''

        with open(filename, 'w') as handle:
            for complex in complexToRoles:
                handle.write('%s\t%s\n' %(complex, self.separator.join(complexToRoles[complex])))
        return

    # The reaction to complexes file contains a mapping of reaction IDs to complex IDs.
    # Each line has these fields:
    #   1. Reaction ID in KBase format (e.g. kb|rxn.5682)
    #   2. List of complex IDs in KBase format (e.g. kb|cpx.1507///kb|cpx.1813)
    #
    # Note that complex IDs must be separated by a string that does not occur in any complex ID.
    
    def readReactionComplexFile(self, filename):
        ''' Read data from a reaction to complexes file.

            @param filename: Path to reaction to complexes file
            @return Dictionary mapping a reaction ID to list of complex IDs
        '''

        rxnToComplexes = dict()
        with open(filename, 'r') as handle:
            for line in handle:
                fields = line.strip('\r\n').split('\t')
                rxn = fields[0]
                complexes = fields[1].split(self.separator)
                if rxn not in rxnToComplexes:
                    rxnToComplexes[rxn] = complexes
                else:
                    rxnToComplexes[rxn] += complexes
        return rxnToComplexes

    def writeReactionComplexFile(self, filename, rxnToComplexes):
        ''' Write data to a reaction to complexes file.

            @param filename: Path to reaction to complexes file
            @param rxnToComplexes: Dictionary mapping a reaction ID to list of complex IDs
            @return Nothing
        '''

        with open(filename, 'w') as handle:
            for rxn in rxnToComplexes:
                handle.write('%s\t%s\n' %(rxn, self.separator.join(rxnToComplexes[rxn])))
        return
    
    # The reaction mapping file contains a mapping of KBase reactions to KEGG reactions.
    # Each line has these fields:
    #   1. Reaction ID in KEGG format (e.g. R03056)
    #   2. Reaction ID in KBase format (e.g. )
    #   3. Name of reaction in CDM

    def readKeggReactionFile(self, filename):
        ''' Read data from a KEGG to KBase reaction mapping file.

            The reaction mapping is a tuple where the first element is the KEGG
            reaction ID, the second element is the KBase reaction ID, and the
            third element is the KBase reaction name.

            @param filename: Path to reaction mapping file
            @return List of tuples as described above
        '''
        keggReactionList = list()
        lastRxn = 'R00000'
        with open(filename, 'r') as handle:
            for line in handle:
                fields = line.strip('\r\n').split('\t')
                if fields[0] == lastRxn:
                    print fields
                else:
                    lastRxn = fields[0]
                keggReactionList.append( [ fields[0], fields[1], fields[2] ] )
        return keggReactionList

    def writeKeggReactionFile(self, filename, keggRxnIdList):
        ''' Write data to KEGG to KBase reactions file.

            The input list is a list of tuples where the first element is the
            KEGG reaction ID, the second element is the KBase reaction ID, and
            the third element is the KBase reaction name.

            @param filename: Path to reaction mapping file
            @param keggRxnIdList: List of tuples as described above
            @return Nothing
        '''

        with open(filename, 'w') as handle:
            for index in range(len(keggRxnIdList)):
                handle.write('%s\t%s\t%s\n' %(keggRxnIdList[index][0], keggRxnIdList[index][1], keggRxnIdList[index][2]))
        return

    def mergeFiles(self, sourceFiles, targetFile):
        ''' Merge a list of source files into a target file.

            @param sourceFiles: List of paths to source files
            @param targetFile: Path to target file
            @return Nothing
        '''

        # Open the target file.
        with open(targetFile, 'w') as target:
            # Write every line from all source files to the target file.
            for srcfile in sourceFiles:
                with open(srcfile, 'r') as source:
                    for line in source:
                        target.write(line)
        return

    def mergeFilesById(self, sourceFiles, targetFile):
        ''' Merge a list of source files into a target file with no duplicate IDs.

            Each source file must have two tab delimited fields where the first
            field is an ID and the second field is a list of values delimited by
            the separator configuration variable.

            @param sourceFiles: List of paths to source files
            @param targetFile: Path to target file
            @return Nothing
        '''

        # Start with an empty output dictionary.
        mergedDict = dict()

        # Read every line from every source file.  If the ID is not a key in the
        # output dictionary, add it to the dictionary.  Otherwise, append the value
        # to the current value using the separator configuration variable.
        for index in range(len(sourceFiles)):
            with open(sourceFiles[index], 'r') as handle:
                for line in handle:
                    fields = line.strip('\r\n').split('\t')
                    key = fields[0]
                    if key not in mergedDict:
                        mergedDict[key] = fields[1]
                    else:
                        value = mergedDict[key]+self.separator+fields[1]
                        mergedDict[key] = value

        # Write the output dictionary to the target file.
        with open(targetFile, 'w') as handle:
            for key in sorted(mergedDict):
                handle.write('%s\t%s\n' %(key, mergedDict[key]))

        return

    def mergeDataFiles(self):
        ''' Merge the data files from the configured sources.

            @note There must be a matching key in the dictionary of files for a
                source as there is in the DataFiles dictionary.
            @return Nothing
        '''

        # The data files are created by merging the data files from each of the
        # configured sources.
        for key in self.DataFiles:
            # Build the list of source files to be merged.
            sourceList = list()
            for src in self.sources:
                if len(self.sources[src]) > 0: # Make sure source is configured
                    sourceList.append(self.sources[src][key])

            # Merge the source files to the target file.
            if key == 'reaction_complex_file' or key == 'complex_role_file':
                self.mergeFilesById(sourceList, self.DataFiles[key])
            else:
                self.mergeFiles(sourceList, self.DataFiles[key])

        return

    def readRolesetProbabilityFile(self, roleset_probability_file):
        ''' Read the roleset probability file.
        
            @param roleset_probability_file Path to roleset probability file
            @return Dictionary mapping query ID to list of tuples with list of roles and probability
        '''
    
        queryToTuplist = dict()
        for line in open(roleset_probability_file, 'r'):
            spl = line.strip('\r\n').split('\t')
            if spl[0] in queryToTuplist:
                queryToTuplist[spl[0]].append( (spl[1], float(spl[2])) )
            else:
                queryToTuplist[spl[0]] = [ (spl[1], float(spl[2])) ]
        return queryToTuplist
    
    # The status file is used to track the status of setting up the static database files when
    # the server starts.  The first line of the file contains the status which is one of
    # these values:
    #   1. 'building' when the pa-gendata command is building the files
    #   2. 'running' when the server initialization is in progress
    #   3. 'ready' when the server initialization is complete or a build is complete
    #   4. 'failed' when there was an error building or loading the files
    #
    # The second line has the timestamp of when the status was last changed.
    
    def readStatusFile(self):
        ''' Read the current status value from the status file.
        
            @return Current status string
        '''
    
        with open(self.StatusFiles['status_file'], 'r') as handle:
            statusLine = handle.readline()
        return statusLine.strip('\r\n')
    
    def writeStatusFile(self, status):
        ''' Write new status value to the status file.
        
            @param status New status value
            @return Nothing
        '''
    
        with open(self.StatusFiles['status_file'], 'w') as handle:
            handle.write('%s\nupdated at %s\n' %(status, time.strftime('%a %b %d %Y %H:%M:%S %Z', time.localtime())))
        return
    
    def checkIfDatabaseFilesExist(self):
        ''' Check for existence of all of the database files.
        
            @raise NotReadyError: A database file does not exist
            @return Nothing
        '''
    
        for path in self.DataFiles.values():
            if not os.path.exists(path):
                raise NotReadyError('Static database file "%s" does not exist' %(path))
        for path in self.SearchFiles.values():
            if not os.path.exists(path):
                raise NotReadyError('Static database file "%s" does not exist' %(path))
        return

    def loadDatabaseFiles(self, mylog):
        ''' Load the static database files from Shock.

            The static database files are stored in the directory specified by the
            data_dir configuration variable.  A file is only downloaded if
            the file is not available on this system or the file has been updated
            in Shock.

            @param mylog Log object for messages
            @return Nothing
            @raise MissingFileError when database file is not found in Shock
        '''
        
        # Get the current info about the static database files from the cache file.
        cacheFilename = self.StatusFiles['cache_file']
        if os.path.exists(cacheFilename):
            fileCache = json.load(open(cacheFilename, 'r'))
        else:
            fileCache = dict()
        
        # Create a shock client.
        shockClient = ShockClient(self.shockURL)

        # See if the static database files on this system are up-to-date with files stored in Shock.
        shockFiles = dict(self.DataFiles.items() + self.SearchFiles.items())
        for key in shockFiles:
            # Get info about the file stored in Shock.
            localPath = shockFiles[key]
            name = os.path.basename(localPath)
            nodelist = shockClient.query_node( { 'lookupname': 'ProbAnnoData/'+name } )
            if len(nodelist) == 0:
                message = 'Database file %s is not available from %s\n' %(name, self.shockURL)
                mylog.log_message(log.ERR, message) # MBM
                raise MissingFileError(message)
            node = nodelist[0]
            
            # Download the file if the checksum does not match or the file is not available on this system.
            download = False
            if key in fileCache:
                if node['file']['checksum']['md5'] != fileCache[key]['file']['checksum']['md5']:
                    download = True
            else:
                download = True
            if os.path.exists(localPath) == False:
                download = True
            if download:
                shockClient.download_to_path(node['id'], localPath)
                fileCache[key] = node
                mylog.log_message(log.INFO, 'Downloaded %s to %s' %(key, localPath))
                
        # Save the updated cache file.
        json.dump(fileCache, open(cacheFilename, 'w'), indent=4)
        return
     
    def storeDatabaseFiles(self, token):
        ''' Store the static database files to Shock.

            @param token: Authorization token for authenticating to shock
            @return Nothing
        '''
        
        # Create a shock client.
        shockClient = ShockClient(self.shockURL, token=token)
        
        # Upload all of the static database files to shock.
        fileCache = dict()
        shockFiles = dict(self.DataFiles.items() + self.SearchFiles.items())
        for key in shockFiles:
            localPath = shockFiles[key]
            name = os.path.basename(localPath)
            if os.path.exists(localPath):
                sys.stderr.write('Saving "%s"...' %(localPath))
                
                # See if the file already exists in Shock.
                query = { 'lookupname': 'ProbAnnoData/'+name }
                nodelist = shockClient.query_node(query)
                
                # Remove all instances of the file in Shock.
                if nodelist != None:
                    for node in nodelist:
                        shockClient.delete_node(node['id'])
     
                # Build the attributes for this file and store as json in a separate file.
                moddate = time.ctime(os.path.getmtime(localPath))           
                attr = { 'lookupname': 'ProbAnnoData/'+name, 'moddate': moddate }
                attrFilename = os.path.join(self.dataFolderPath, name+'.attr')
                attrFid = open(attrFilename, 'w')
                json.dump(attr, attrFid, indent=4)
                attrFid.close()
                
                # Upload the file to Shock.
                metadata = shockClient.create_node(localPath, attrFilename)
                fileCache[key] = metadata
                os.remove(attrFilename)
                
                # Remove the list of users from the read ACL to give the file public read permission.
                # Note this needs to change for Shock version 0.9.5 but not sure how to set public ACLs.
                readacl = shockClient.get_acl(metadata['id'])
                shockClient.delete_acl(metadata['id'], 'read', readacl['read'][0])
                sys.stderr.write('done\n')
                
            else:
                sys.stderr.write('Could not find "%s" so it was not saved\n' %(localPath))
                
        # Save the metadata on all of the database files.
        cacheFilename = os.path.join(self.dataFolderPath, StatusFiles['cache_file'])
        json.dump(fileCache, open(cacheFilename, 'w'), indent=4)

        return

    def getDatabaseFiles(self, mylog, testDataPath):
        ''' Get the static database files.

            The static database files come from one of three places: (1) Shock,
            (2) a local data directory, (3) a test data directory.  If the files
            are not found in Shock, the local data directory is searched.  If the
            file are not found in the local data directory, the test data directory
            is used.

            @param mylog: Log object for messages
            @param testDataPath: Path to directory with test database files
            @return Current value of load data option which indicates which of the
                three places is being used for the static database files
        '''

        # Update the status file to indicate that the static database files are being updated.
        self.writeStatusFile('running')
        status = 'failed'

        # Get the static database files from Shock (only missing or changed files are downloaded).
        if self.loadDataOption == 'shock':
            try:
                self.loadDatabaseFiles(mylog)
                status = 'ready'
                mylog.log_message(log.INFO, 'All static database files loaded from Shock to %s' %(self.dataFolderPath))
            except:
                traceback.print_exc(file=sys.stderr)
                mylog.log_message(log.NOTICE, 'Failed to load static database files from Shock. Checking current files...')
                self.loadDataOption = 'preload'

        # Get the static database files from the data directory specified in the configuration.
        if self.loadDataOption == 'preload':
            try:
                self.checkIfDatabaseFilesExist()
                status = 'ready'
                mylog.log_message(log.INFO, 'All static database files are available in %s' %(self.dataFolderPath))
            except:
                # There is a problem with at least one of the static database files so switch
                # to the test data.
                status = 'ready'
                self.loadDataOption = 'test'
                self.dataFolderPath = testDataPath
                traceback.print_exc(file=sys.stderr)
                mylog.log_message(log.NOTICE, 'Static database files are missing. Switched to test database files in %s' %(testDataPath))

        # Update the status file to indicate that the static database files updating is done.
        self.writeStatusFile(status)
        return self.loadDataOption
