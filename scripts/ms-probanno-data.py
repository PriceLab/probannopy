#! /usr/bin/env python

import argparse
import os
import sys
import traceback
from ConfigParser import ConfigParser
from biop3.ProbModelSEED.ProbAnnotationParser import ProbAnnotationParser
from biokbase import log

desc1 = '''
NAME
      ms-probanno-data -- manage data files for probabilistic annotation algorithm

SYNOPSIS
'''

desc2 = '''
DESCRIPTION
      Manage the data files required by the probabilistic annotation algorithm.
      
      The action argument specifies the action to perform. The following actions
      are supported: (1) "load" to load the data files from Shock, (2) "store"
      to store the data files to Shock, or (3) "builddb" to build a search
      database for the configured search program.
       
      The --token optional argument specifies the authentication token for the
      user and is required when using the "store" action.
      
      The Shock service URL, data directory, and search program are obtained
      from the deployed configuration file.
'''

desc3 = '''
EXAMPLES
      Load the static data files from Shock:
      > ms-probanno-data load

AUTHORS
      Mike Mundy 
'''

if __name__ == '__main__':
    # Parse options.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='ms-probanno', epilog=desc3)
    parser.add_argument('action', help='action to perform (load, store, builddb)', action='store', default=None)
    parser.add_argument('--token', help='token for user', action='store', dest='token', default=None)
    usage = parser.format_usage()
    parser.description = desc1 + '      ' + usage + desc2
    parser.usage = argparse.SUPPRESS
    args = parser.parse_args()
    
    # Get the configuration variables.
    serviceName = os.environ.get('KB_SERVICE_NAME', 'ProbModelSEED')
    cfg = ConfigParser()
    cfg.read(os.path.join(os.environ.get('KB_TOP'), 'deployment.cfg'))
    config = dict()
    for nameval in cfg.items(serviceName):
        config[nameval[0]] = nameval[1]
        
    # Create a logger.
    logger = log.log(serviceName, ip_address=True, authuser=True, module=True, method=True,
        call_id=True, logfile=config['mlog_log_file'])
    logger.set_log_level(int(config['mlog_log_level']))

    # Create a ProbAnnotationParser object for working with the static data files.
    dataParser = ProbAnnotationParser(config)
    
    # Load the static data files from Shock.
    if args.action == 'load':
        print 'Started loading data files from Shock ...'
        dataParser.writeStatusFile('running')
        try:
            dataParser.loadDatabaseFiles(logger)
            dataParser.writeStatusFile('ready')
        except:
            dataParser.writeStatusFile('failed')
            print 'Failed to load static data files from Shock'
            traceback.print_exc(file=sys.stderr)
            exit(1)
    
    # Store the static data files to Shock.
    elif args.action == 'store':
        print 'Started storing data files to Shock ...'
        try:
            dataParser.storeDatabaseFiles(args.token)
        except:
            print 'Failed to store static data files to Shock'
            traceback.print_exc(file=sys.stderr)
            exit(1)            
    
    # Build a search database for the configured search program.
    elif args.action == 'builddb':
        print 'Started building search database ...'
        try:
            dataParser.buildSearchDatabase()
        except MakeblastdbError as e:
            print 'Failed to build a search database: '+e.message
        
    else:
        print 'Action '+args.action+' is not supported'
        exit(1)
        
    exit(0)
