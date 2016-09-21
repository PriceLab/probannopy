#! /usr/bin/env python
# BASE_PATH is the absolute path of .. relative to this script location
import sys
import os
BASE_PATH = reduce (lambda l,r: l + os.path.sep + r,
        os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep )[:-1] )
sys.path.append( os.path.join( BASE_PATH, "lib" ) )
sys.path.append( os.path.join( BASE_PATH, "lib/python2.6/site-packages" ) )
sys.path.append( os.path.join( BASE_PATH, "../lib/python2.6/site-packages" ) )
import argparse
import json
import traceback
import re
import urllib2
import random
import datetime
import wget

from ProbAnnotationWorker import ProbAnnotationWorker

desc1 = '''
NAME
      ms-probanno-standalone -- run probabilistic annotation algorithm for a genome

SYNOPSIS
'''

desc2 = '''
DESCRIPTION
      Run the probabilistic algorithm for a genome.  The proteome argument is
      a protein fasta file or Uniprot ID.  The templatefile argument is
      a json file containing the template model used to reconstruct a model for the
      organism.  The rxnprobsfile argument is the file path to where the output
      rxnprobs object will be stored.
'''

desc3 = '''
EXAMPLES
      Run probabilistic annotation for Bacillus subtilis genome:
      > ms-probanno-standalone genomes/1415167.3.PATRIC.faa
          templates/GramNegative.json 1415167.3.probanno.out

AUTHORS
      Mike Mundy, Terry Farrah
'''


def writeRxnprobs(rxnProbs, filename, templatefile, genome_id, organism):
    """ Write a tab-delimited file of reaction probabilities
        reactionProbs is a list of
        rxn (string), maxProb, TYPE (string), complexString, GPR (string)
    """
    try:
	f = open(filename, 'w')
    except:
	sys.exit(2) #exit code 2: cannot open output file for writing
    f.write("# ProbAnno run " + str(datetime.datetime.now()) + "\n")
    f.write("# genomeID:" + genome_id + " template:" + os.path.splitext(os.path.basename(templatefile))[0] + " organism:" + organism + "\n")
    if len(rxnProbs) == 0:
	sys.exit(3)  # exit code 3: empty results
    for index in range(len(rxnProbs)):
	prob = rxnProbs[index]
	f.write('{0}0\t{1:1.4f}\t{2}\t{3}\t{4}\n'.format(prob[0], prob[1], prob[2], prob[3], prob[4]))
    f.close()
    return


if __name__ == '__main__':
    # Parse options.
#f = open("/local/local_webservices/probanno/test.log", 'w')
#f.write( str(datetime.datetime.now()) + "\n")
#f.close()
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     prog='ms-probanno-standalone', epilog=desc3)

    # Get proteome -- either a fasta file or a Uniprot ID
    parser.add_argument('proteome',
                        help='fasta file of protein sequences OR taxonomy identifier for proteome (e.g. 224308 for Bacillus subtilis)',
                        action='store', default=None)

    # Get template file, previously created using Build_Model_Template.py from ModelSEEDDatabase
    parser.add_argument('templatefile', help='json template file generated using Build_Model_Template.py',
                        action='store', default=None)

    # Get output filename
    parser.add_argument('rxnprobsfile', help='output filename for rxnprobs', action='store', default=None)

    # Get optional genome_id
    cl_genome_id = ''
    parser.add_argument('--genome_id', dest='genome_id',  help='a string denoting this genome, often an organism name')

    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()
    random.seed()
    fastaFile = args.proteome
    genome_id = os.path.splitext(os.path.basename(fastaFile))[0]
    if args.genome_id:
            genome_id = args.genome_id
	    print "args.genome_id: " + args.genome_id
    

    # Create a worker for running the algorithm.
    worker = ProbAnnotationWorker(genome_id)
        
    # taxon id is one to seven digits
    taxid_pattern = re.compile('^\d{1,7}$')
    if taxid_pattern.match(args.proteome):  # fetch file from Uniprot

        url = "http://www.uniprot.org/uniprot/?format=fasta&query=organism:" + args.proteome
        fastaFile = worker.workFolder + "/" + args.proteome + "_" + str(random.randint(100000, 999999)) + ".fasta#"
        attempts = 0
        while attempts < 3:
            try:
		print "Fetching " +  url
                response = urllib2.urlopen(url, timeout = 5)
                content = response.read()
                f = open( fastaFile, 'w' )
                f.write( content )
                f.close()
                break
            except urllib2.URLError as e:
                attempts += 1
                print type(e)
		sys.exit(7)     # exit code 7: could not fetch fasta file
	if os.path.getsize(fastaFile) < 10:
		sys.exit(8)	# exit code 8: empty fasta file
    elif args.proteome.startswith("PATRIC:"):
        patric_genome_id = args.proteome[7:]
        url = "ftp://ftp.patricbrc.org/patric2/patric3/genomes/" + patric_genome_id + "/" + patric_genome_id + ".PATRIC.faa"
        fastaFile = wget.download(url)
        os.rename(fastaFile, "genomes/" + fastaFile)
        fastaFile = "genomes/" + fastaFile

   # Get organism name from fasta file
    f = open( fastaFile, 'r' )
    line = f.readline()
    org_pattern = re.compile('.*OS=(.+) GN=')
    match = org_pattern.match(line)
    args.organism = 'unspecified organism'
    if match:
	    args.organism = match.group(1)
    f.close()

    # Create a dictionary from the json template file
    json_data = open(args.templatefile).read()
    template = json.loads(json_data)

    # Build a dictionary to look up roles in the template by ID.
    roles = dict()
    for index in range(len(template['roles'])):        
        roles[template['roles'][index]['id']] = index

    # Create a dictionary to map a complex to a list of roles as defined in the template.
    complexesToRoles = dict()
    for index in range(len(template['complexes'])):
        complexId = template['complexes'][index]['id']
        if len(template['complexes'][index]['complexroles']) > 0:
            complexesToRoles[complexId] = list()
            for crindex in range(len(template['complexes'][index]['complexroles'])):
                # A complex has a list of complexroles and each complexrole has a reference
                # to a role and each role has a name. Role ID is last element in reference.
                roleId = template['complexes'][index]['complexroles'][crindex]['templaterole_ref'].split('/')[-1]
                # Mike Mundy used this statement, but it results in all probs being zero.
                #complexesToRoles[complexId].append(roleId)
		complexesToRoles[complexId].append(template['roles'][roles[roleId]]['name']) 

    # Create a dictionary to map a reaction to a list of complexes as defined in the template.
    reactionsToComplexes = dict()
    for index in range(len(template['reactions'])):
        reactionId = template['reactions'][index]['id']
        if len(template['reactions'][index]['templatecomplex_refs']) > 0:
            reactionsToComplexes[reactionId] = list()
            for complexRef in template['reactions'][index]['templatecomplex_refs']:
                # Complex ID is last element in reference.
                reactionsToComplexes[reactionId].append(complexRef.split('/')[-1])

    # Run the probabilistic annotation algorithm.
    try:

        # Run blast using the fasta file.
        blastResultFile = worker.runBlast(fastaFile)
        if taxid_pattern.match(args.proteome):
            os.remove(fastaFile)
        
        # Calculate roleset probabilities.
        rolestringTuples = worker.rolesetProbabilitiesMarble(blastResultFile)
        
        # Calculate per-gene role probabilities.
        roleProbs = worker.rolesetProbabilitiesToRoleProbabilities(rolestringTuples)
        
        # Calculate whole cell role probabilities.
        totalRoleProbs = worker.totalRoleProbabilities(roleProbs)

        # Calculate complex probabilities.
        complexProbs = worker.complexProbabilities(totalRoleProbs, complexesToRequiredRoles=complexesToRoles)
 
        # Calculate reaction probabilities.
        reactionProbs = worker.reactionProbabilities(complexProbs, rxnsToComplexes=reactionsToComplexes)

        # Cleanup work directory.
        worker.cleanup()
        
    except Exception as e:
        worker.cleanup()
        sys.stderr.write('Failed to run probabilistic annotation algorithm: %s\n' % e.message)
        tb = traceback.format_exc()
        sys.stderr.write(tb)
        sys.exit(1)   # exit code 1: could not run ProbAnno algorithm

    # Create output file with details on reaction likelihoods.
    # reactionProbs is a list of rxn (string), maxProb, TYPE (string), complexString, GPR (string)
    writeRxnprobs(reactionProbs, args.rxnprobsfile, args.templatefile, genome_id, args.organism)  # first arg is a list

    sys.exit(0)
