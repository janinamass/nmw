====
nmw
=====

*nmw* is an implementation of the Needleman-Wunsch pairwise global alignment algorithm in python.
All sequences in file A will be compared to all sequences in file B.

Usage ::
	
	nmw.py -m MAT -f FILE_A -g FILE_B -o OUTFILE


  general options:
    -m, --matrix=MAT    (see matrix options)                 
    -f, --fasta_a=FILE    multifasta file a
    -g, --fasta_b=FILE    multifasta file b
    -o, --outfile=FILE    store results in FILE
    -h, --help            prints this
    -n, --num_cores=NUM   number of processors to use [default 1]

  alignment options:
     -O, --gap_open=FLOAT           gap opening cost [default 10]
     -E, --gap_extend=FLOAT         gap extension cost [default 0.5]

  fasta options:
    -d, --delim=STRING               split fasta headers at STRING
    -a, --asID=INT                   use INTth part of fasta header as transcript-ID
                                     (default:0)
  matrix options:
    MAT can be [BLOSUM62, ... todo]

Installation
============
::
    pip install nmw.py


Requirements
============

Python modules
~~~~~~~~~~~~~~~
* pickle

::   
    pip install pickle

