**New!**: the current implementation of `probannopy` is compatible with Python 3.7+. Legacy support for python 2.7.X is maintained on the `python2` branch

Mike Mundy, Matt Benedict, Terry Farrah, and Brendan King March 2017
Probabilistic annotation for metabolic modeling.
Code is derived from git repository ProbModelSEED.

Probanno can be used in a few ways:
 - As a python package. If you want to write a python tool that makes programatic use of probanno, this is the simplest way.
 - As a stand-alone script. See [Probanno-Standalone](https://github.com/tfarrah/ProbAnno-Standalone)
 - From a hosted web service. No installation needed! (Coming Soon)
 
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

* Check that python3 and the required modules are installed. This can be done by navigating to the probanno directory and installing all requirements with pip: `pip install -r requirements.txt`

* Install usearch binary
  (http://www.drive5.com/usearch/manual/install.html)

* Select a directory with at least 1.4G available space (use probanno/data if
  enough space there is the default) and download two data files to that location:
  
    - [PROTEIN.udb](https://drive.google.com/file/d/0B3QgVGEsPx9kS3Y3WkNuSi02ams/view?usp=sharing)
    - [OTU_FID_ROLE](https://drive.google.com/file/d/0B3QgVGEsPx9keW9PYUhDTFFNWWc/view?usp=sharing)
    - If you use another folder than `probanno/data`, then edit `probanno/deploy.cfg` to include the directive: `data_dir=/path/to/folder/containing/data/files/`

* Edit the config file (`probanno/deploy.cfg`) to set the `search_program_path` variable to where you installed
  usearch binary

* Install dependencies
    - (from the directory this repository is downloaded to): `pip install -r probanno/requirements.txt`

* Install probanno with `pip install -e`
    - (from the directory this repository is downloaded to): `pip install -e ./probanno`
    
* To use, `import probanno` in a python shell or script and access functions!

Example invocation:
```python
import probanno
probanno.generate_reaction_probabilities('my_fasta_file.fasta', 'templates/GramNegative.json', genome_id='my_genome')
```

First argument is proteome fasta file for organism of interest. You can retrieve then by PATRIC genome ID using `probanno.get_fasta_by_id` (search for genomes from the [PATRIC home page](https://www.patricbrc.org/portal/portal/patric/Home)) OR Uniprot ID for proteome.
Second argument is appropriate template file of several available. (Included in `/templates`, e.g. `/templates/GramNegative.json`)

