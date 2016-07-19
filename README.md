Mike Mundy, Matt Benedict, and Terry Farrah June 2016
Probabilistic annotation for metabolic modeling.
Code is derived from git repository ProbModelSEED.
Dependencies on servers have been removed.

Originally described in Benedict, et al., Likelihood-Based Gene Annotations for
Gap Filling and Quality Assessment in Genome-Scale Metabolic Models, PLos
Comput Biol October 16, 2014, PMID: 25329157 
Please cite this reference in any work that makes use of this software.

Preparation (commands given are for bash shell):

- Install usearch to a place in $PATH (http://www.drive5.com/usearch/manual/install.html)

- Define the environment variable KB_TOP to the parent of the current directory

  export KB_TOP=..

- Add $KB_TOP/ProbAnno-Standalone/lib to $PYTHONPATH

  export PYTHONPATH=$PYTHONPATH:$KB_TOP/ProbAnno-Standalone/lib

- Select a directory with at least 1.4G available space, download two data
  files to that location, then symlink back to ProbAnno-Standalone directory.

  DATADIR=/foo/bar/baz
  cd $DATADIR
  wget https://www.dropbox.com/s/lucq1p7zd9mmf1j/OTU_FID_ROLE?dl=0 -O OTU_FID_ROLE
  wget https://www.dropbox.com/s/bssrfllefzvhzvu/PROTEIN.udb?dl=0 -O PROTEIN.udb
  cd KB_TOP/ProbAnno-Standalone/data
  ln -s $DATADIR/OTU_FID_ROLE .
  ln -s $DATADIR/PROTEIN.udb .


Example invocation:
scripts/ms-probanno-standalone.py genomes/1415167.3.PATRIC.faa templates/1415167.3.PATRIC.faa 1415167.3.probanno.out
