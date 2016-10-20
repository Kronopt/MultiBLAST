# MultiBLAST
Sequentially Blast every fasta file (.faa or .fas) on a folder against a local database

#### Dependencies
* Python 2.7
* Blast+ executables (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* Linux (tested and developed on Ubuntu 14.0.4, but will most likely work on WINDOWS and MAC too)

#### How to run
* Organize your files as follows on a folder:
  * MultiBLAST.py
  * database files (obtained via makeblast, with `makeblastdb -in [yourGenes].pfasta -dbtype prot`)
  * folder with fasta files (.faa or .fas) to blast (must be the only existing folder in the directory)
* Run MultiBLAST.py from the command-line as follows `python MultiBLAST.py [-t [#]]`
* For more help, type `python MultiBLAST.py -h`

#### Note
I developed this script to help me with a project I was involved in. As much as I tried to make it work on many different scenarios, it is still probably biased to work with the files I was using. Nevertheless, it should work fine for protein blasts. If you want to call more blast+ parameters, edit the variable `standard_blast` inside the `loop_blast` function.
