Mike Mundy, Matt Benedict, and Terry Farrah June 2016
Probabilistic annotation for metabolic modeling.
Code is derived from git repository ProbModelSEED
(https://github.com/ModelSEED/ProbModelSEED), with
dependencies on servers removed.

SOFTWARE AND DATA ATTRIBUTION
-----------------------------

This product includes software developed by and/or derived from the SEED Project
(http://www.theseed.org/) to which the U.S. Government retains certain rights.
Please see the original software license
(https://github.com/ModelSEED/ProbModelSEED/blob/master/LICENSE.md) for details.

The source for the two data files is the KBase Central Data Model and KEGG.
Details on the KBase data policy and sources are listed on
http://kbase.us/data-policy-and-sources/. In particular, ProbAnno-Standalone users
must abide by the following KBase Data Sharing Policy:

KBase conforms to the Information and Data Sharing Policy of the Genomic
Science Program of the Office of Biological and Environmental Research within
the Office of Science. This requires that all publishable data, metadata, and
software resulting from research funded by the Genomic Science program must
conform to community-recognized standard formats when they exist; be clearly
attributable; and be deposited within a community-recognized public database(s)
appropriate for the research.

The KEGG data were obtained using the KEGG REST API. KEGG is described here:
* Kanehisa, et al., KEGG as a reference resource for gene and protein annotation.
Nucleic Acids Res. 44 (2016). PMID 26476454
* Kanehisa and Goto, KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids
Res. 28, 27-30 (2000). PMID: 10592173 

Computational method originally described in:
* Benedict, et al., Likelihood-Based Gene Annotations for Gap Filling and Quality
Assessment in Genome-Scale Metabolic Models. PLos Comput Biol October 16, 2014,
PMID: 25329157 

Please cite all three references above in any work that makes use of this software.


INSTALLATION and USE
--------------------

* Check that python and the required modules are installed.
  Required modules are  argparse, json, sys, os, time, traceback, requests, re,
  wget, logging, subprocess, shutil, math, tempfile, datetime, urllib2,
  ConfigParser

* Install usearch
  (http://www.drive5.com/usearch/manual/install.html)

* Create a environment variable for the directory that this README.md file is in

```
export PADIR=$(pwd)   # bash shell; modify if using csh/ksh etc.
```

* Select a directory with at least 1.4G available space (use $PADIR/data if
  enough space there) and download two data files to that location.

```
  bash
  DATADIR=/foo/bar/baz  # fill in appropriate path. Again, this is for bash shell.
  cd $DATADIR
  wget -O OTU_FID_ROLE https://www.dropbox.com/s/lucq1p7zd9mmf1j/OTU_FID_ROLE?dl=0
  wget -O PROTEIN.udb https://www.dropbox.com/s/bssrfllefzvhzvu/PROTEIN.udb?dl=0
```
* Edit the config file to set the `search_program_path` variable to where you installed
  usearch and the `data_dir` variable to $DATADIR:

```
  bash
  vim $PADIR/deploy.cfg
```



Example invocation:


First argument is proteome fasta file for organism of interest OR PATRIC genome ID (search for genomes from the [PATRIC home page](https://www.patricbrc.org/portal/portal/patric/Home)) OR Uniprot ID for proteome.
Second argument is appropriate template file of several available.
Third argument is output filename.

```bash
cd $PADIR
scripts/ms-probanno-standalone.py genomes/1415167.3.PATRIC.faa templates/GramNegative.json 1415167.3.probanno.out # silent; takes a few minutes.
scripts/ms-probanno-standalone.py PATRIC:226186.12 templates/GramNegative.json 226186.12.probanno.out # silent; takes a few minutes.
```
