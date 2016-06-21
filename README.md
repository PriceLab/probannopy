Terry Farrah June 2016
Probabilistic annotation for metabolic modeling.
Code is derived from git repository ProbModelSEED.
Dependencies on servers have been removed.

Preparation:
- Install usearch to a place in $PATH (http://www.drive5.com/usearch/manual/install.html)
- Define the environment variable KB_TOP to the parent of the current directory
export KB_TOP=..
- Add $KB_TOP/ProbAnno-Standalone/lib to $PYTHONPATH
export PYTHONPATH=$PYTHONPATH:$KB_TOP/ProbAnno-Standalone/lib

Example invocation:
scripts/ms-probanno.py genomes/1415167.3.PATRIC.faa templates/1415167.3.PATRIC.faa output
