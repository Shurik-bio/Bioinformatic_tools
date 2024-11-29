# Bioinformatics tools

`Bioinformatics_tools` utility allows analyze nucleotide and aminoacids sequences with help `run_dna_rna_tools.py` as well as `filter_fastq.py` filters raw FASTA reads according to their quality and set parameters.
With using functions from module `dna_rna_tools.py` you can return transcribe, reverse, complement or reverse complement nucleotide sequences. 

The `Bioinformatics_tools` package has the following structure:

 ```python
    -/
     |- README.md
     |- main_script.py 
     |- bio_files_processor.py 
     |- modules/
           |- dna_rna_tools.py
           |- filter_fastq.py
           |- ...

```

Example of using `run_dna_rna_tools`, which is acceptable for one or more input sequences:
```
run_dna_rna_tools('CGT', 'tgA', 'transcribe') # ['CGU', tgU] 
run_dna_rna_tools('GTTAC', 'reverse') # 'CATTG'
run_dna_rna_tools('TGaT', 'complement') # 'ACtA'
run_dna_rna_tools('ACtA', 'reverse_complement') # 'TaGT'

 The complement of a sequence is given by replacing each base with its complementary base, as given by the dictionary:
    {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'U': 'A', 'u': 'a', 'a': 't', 'c': 'g', 't': 'a', 'g': 'c'}
```

`Filter_fastq` module filters a FASTQ file based on GC%, length and quality thresholds and writes the filtered sequences to a new fastq file. It contains such functions as:

- `gc_bounds` takes a list of integers representing GC% borders and filters a fastq file based on those borders. If the list contains one element, it is considered as a high GC% border and 
    all sequences with GC% below it are filtered out. If the list contains two elements, the function filters out all sequences with GC% below the first element and above the second one.

-  `length_bounds` takes a list of integers representing length borders and filters a fastq file based on those borders. If the list contains one element, it is considered as a low length 
     border and all sequences with length below it are filtered out. If the list contains two elements, the function filters out all sequences with length below the first element and above the second one.

-  `quality_threshold` takes an integer representing a quality threshold and filters a fastq file based on that threshold. The function calculates the quality of each read in the file, and if 
     the quality is above the given threshold, the read is included in the filtered dictionary.


 `Bio_files_processor.py` utility converts a multi-line FASTA file to a one-line fasta file and parses a BLAST output file and writes the results to an output file.

Built-in function `Converts a multi-line fasta file` accepts the FASTA file with sequences of DNA, RNA or aminoacids, which are splitted on a several lines. It creates a new file in FASTA format, where each sequence is fitted in a single line.  

`Parse_blast_output` utility assepts the BLAST result as input in txt format and save the list of the most significant nucleotide or aminoacids sequences in alphabetical order in output file.











