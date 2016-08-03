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
import time
import traceback
import requests
import re
import urllib2
import random


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

def writeRxnprobs(rxnProbs, file):
    ''' Write a tab-delimited file of reaction probabilities
        reactionProbs is a list of
	rxn (string), maxProb, TYPE (string), complexString, GPR (string)
    '''
    f = open(file, 'w')
    while rxnProbs:
        item = rxnProbs.pop(0)
        s = str(item.pop(0))
	for i in range (1, 5):
	    element = str(item.pop(0))
	    s = "\t".join([s,element])
        s = s + "\n"
	f.write(s)


if __name__ == '__main__':
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='ms-probanno', epilog=desc3)

    # Get proteome -- either a fasta file or a Uniprot ID
    parser.add_argument('proteome', help='fasta file of protein sequences OR uniprot identifier for proteome (e.g. UP000001570)', action='store', default=None)

    # Get template file, previously created using Build_Model_Template.py from ModelSEEDDatabase
    parser.add_argument('templatefile', help='json template file generated using Build_Model_Template.py', action='store', default=None)

    # Get output filename
    parser.add_argument('rxnprobsfile', help='output filename for rxnprobs', action='store', default=None)

    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()
    random.seed()
    fastaFile = args.proteome
    genome_id = os.path.splitext(os.path.basename(fastaFile))[0]

    # Create a worker for running the algorithm.
    worker = ProbAnnotationWorker(genome_id)
        
    p = re.compile('^UP0000\d\d\d\d\d$')
    if p.match(args.proteome):  # fetch file from Uniprot
        url = "http://www.uniprot.org/proteomes/" + args.proteome + ".fasta"
        fastaFile = worker.workFolder + "/" + args.proteome + str(random.randint(100000, 999999)) + ".fasta"
	attempts = 0
	while attempts < 3:
	    try:
		response = urllib2.urlopen(url, timeout = 5)
		content = response.read()
		f = open( fastaFile, 'w' )
		f.write( content )
		f.close()
		break
	    except urllib2.URLError as e:
		attempts += 1
		print type(e)

    # Create a dictionary from the json template file
    json_data=open(args.templatefile).read()
    template=json.loads(json_data)

    # Build a dictionary to look up roles in the template by ID.
    roles = dict()
    for index in range(len(template['roles'])):        
        roles[template['roles'][index]['id']] = index

    # Create a dictionary to map a complex to a list of roles as defined in the template.
    complexesToRoles = dict()
    for index in range(len(template['complexes'])):
        complexId = template['complexes'][index]['id']
        complexesToRoles[complexId] = list()
        for crindex in range(len(template['complexes'][index]['complexroles'])):
            # A complex has a list of complexroles and each complexrole has a reference
            # to a role and each role has a name. Role ID is last element in reference.
            roleId = template['complexes'][index]['complexroles'][crindex]['templaterole_ref'].split('/')[-1]
            complexesToRoles[complexId].append(template['roles'][roles[roleId]]['name'])        

    # Create a dictionary to map a reaction to a list of complexes as defined in the template.
    reactionsToComplexes = dict()
    for index in range(len(template['reactions'])):
        reactionId = template['reactions'][index]['id']
        reactionsToComplexes[reactionId] = list()
        for complexRef in template['reactions'][index]['templatecomplex_refs']:
            # Complex ID is last element in reference.
            reactionsToComplexes[reactionId].append(complexRef.split('/')[-1])

    # Run the probabilistic annotation algorithm.
    try:

        # Run blast using the fasta file.
        blastResultFile = worker.runBlast(fastaFile)
        if p.match(args.proteome):
            os.remove(fastaFile)
        
        # Calculate roleset probabilities.
        rolestringTuples = worker.rolesetProbabilitiesMarble(blastResultFile)
        
        # Calculate per-gene role probabilities.
        roleProbs = worker.rolesetProbabilitiesToRoleProbabilities(rolestringTuples)
        
        # Calculate whole cell role probabilities.
        totalRoleProbs = worker.totalRoleProbabilities(roleProbs)

        # Calculate complex probabilities.
        complexProbs = worker.complexProbabilities(totalRoleProbs, complexesToRequiredRoles = complexesToRoles)
 
        # Calculate reaction probabilities.
        reactionProbs = worker.reactionProbabilities(complexProbs, rxnsToComplexes = reactionsToComplexes) 
        # Cleanup work directory.
        worker.cleanup()
        
    except Exception as e:
        worker.cleanup()
        sys.stderr.write('Failed to run probabilistic annotation algorithm: %s\n' %(e.message))
        tb = traceback.format_exc()
        sys.stderr.write(tb)
        exit(1)

    # Create the rxnprobs object in the workspace.
    data = dict()
    data['reaction_probabilities'] = reactionProbs
    # reactionProbs is a list of rxn (string), maxProb, TYPE (string), complexString, GPR (string)
    writeRxnprobs(reactionProbs, args.rxnprobsfile)  # first arg is a list

    exit(0)
