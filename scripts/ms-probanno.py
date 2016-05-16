#! /usr/bin/env python

import argparse
import json
import sys
import os
import time
import traceback
import requests
from biop3.ProbModelSEED.ProbAnnotationWorker import ProbAnnotationWorker
from biop3.Workspace.WorkspaceClient import Workspace, ServerError as WorkspaceServerError, _read_inifile

desc1 = '''
NAME
      ms-probanno -- run probabilistic annotation algorithm for a genome

SYNOPSIS
'''

desc2 = '''
DESCRIPTION
      Run the probabilistic algorithm for a genome.  The genomeref argument is
      the reference to the input genome object.  The templateref argument is
      the reference to the template model used to reconstruct a model for the
      organism.  The rxnprobsref argument is the reference to where the output
      rxnprobs object is stored.
      
      The --ws-url optional argument specifies the url of the workspace service
      endpoint.  The --token optional argument specifies the authentication
      token for the user.
'''

desc3 = '''
EXAMPLES
      Run probabilistic annotation for Bacillus subtilis genome:
      > ms-probanno /mmundy/home/models/.224308.49_model/224308.49.genome
          /chenry/public/modelsupport/templates/GramPositive.modeltemplate
          /mmundy/home/models/.224308.49_model/224308.49.rxnprobs

AUTHORS
      Mike Mundy 
'''

def getObject(wsClient, reference, token):
    ''' Get an object from the workspace.
    
        @param wsClient: Workspace client object
        @param reference: Reference to workspace object
        @return Object data in JSON format
    '''

    requests.packages.urllib3.disable_warnings()
    retryCount = 3
    while retryCount > 0:
        try:
            # The get() method returns an array of tuples where the first element is
            # the object's metadata (which is valid json) and the second element is
            # the object's data (which is not valid json).
            object = wsClient.get({ 'objects': [ reference ] })
            if len(object[0][0][11]) > 0:
                response = requests.get(object[0][0][11]+'?download', headers={ 'AUTHORIZATION': 'OAuth '+token })
                if response.status_code != requests.codes.OK:
                    response.raise_for_status()
                object[0][1] = response.text
            return json.loads(object[0][1])
    
        except WorkspaceServerError as e:
            # When there is a network glitch, wait a second and try again.
            if 'HTTP status: 503 Service Unavailable' in e.message or 'HTTP status: 502 Bad Gateway' in e.message:
                retryCount -= 1
                time.sleep(1)
            else:                
                sys.stderr.write('Failed to get object using reference %s\n' %(reference))
                tb = traceback.format_exc()
                sys.stderr.write(tb)
                exit(1)

    sys.stderr.write('Failed to get object using reference %s because of network problems\n' %(reference))
    exit(1)

def putObject(wsClient, reference, type, data):
    ''' Put an object to the workspace.
    
        @param wsClient: Workspace client object
        @param reference: Reference to workspace object
        @param type: Type of object
        @param data: Object data in JSON format
        @return Nothing
    '''

    retryCount = 3
    while retryCount > 0:
        try:
            wsClient.create({ 'objects': [ [ reference, type, { }, data ] ], 'overwrite': 1 })
            return
    
        except WorkspaceServerError as e:
            # When there is a network glitch, wait a second and try again.
            if 'HTTP status: 503 Service Unavailable' in e.message or 'HTTP status: 502 Bad Gateway' in e.message:
                retryCount -= 1
                time.sleep(1)
            else:                
                sys.stderr.write('Failed to create object using reference %s\n' %(reference))
                tb = traceback.format_exc()
                sys.stderr.write(tb)
                exit(1)
                
    sys.stderr.write('Failed to create object using reference %s because of network problems\n' %(reference))
    exit(1)

if __name__ == '__main__':
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='ms-probanno', epilog=desc3)
    parser.add_argument('genomeref', help='reference to genome object', action='store', default=None)
    parser.add_argument('templateref', help='reference to template model object', action='store', default=None)
    parser.add_argument('rxnprobsref', help='reference to rxnprobs object', action='store', default=None)
    parser.add_argument('--ws-url', help='url of workspace service endpoint', action='store', dest='wsURL', default='https://p3.theseed.org/services/Workspace')
    parser.add_argument('--token', help='token for user', action='store', dest='token', default=None)
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()
    
    # Get the token from the config file if one is not provided.
    if args.token is None:
        authdata = _read_inifile()
        args.token = authdata['token']
     
    # Workaround for extraneous delimiter tacked on the end of references.
    if args.genomeref[-2:] == '||':
        args.genomeref = args.genomeref[:-2]
    if args.templateref[-2:] == '||':
        args.templateref = args.templateref[:-2]
    
    # Get the genome object from the workspace (for the features).
    wsClient = Workspace(url=args.wsURL, token=args.token)
    genome = getObject(wsClient, args.genomeref, args.token)

    # Get the template object from the workspace (for the complexes and roles).
    template = getObject(wsClient, args.templateref, args.token)

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

    # Create a worker for running the algorithm.
    worker = ProbAnnotationWorker(genome['id'])
        
    # Run the probabilistic annotation algorithm.
    try:
        # Convert the features in the genome object to a fasta file.
        fastaFile = worker.genomeToFasta(genome['features'])
        
        # Run blast using the fasta file.
        blastResultFile = worker.runBlast(fastaFile)
        
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
    putObject(wsClient, args.rxnprobsref, 'rxnprobs', data)

    exit(0)
