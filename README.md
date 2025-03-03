# Bioinformatics tools

`Bioinformatics_tools` utility allows analyze nucleotide and aminoacids sequences with help `NucleicAcidSequence` as well as `FastQFilter` filters raw FASTA reads according to their quality and set parameters.
With using functions from module `NucleicAcidSequence` you can return transcribe, reverse, complement or reverse complement nucleotide sequences and check if input sequences related to correct nucleotide or amino acids. 

The `Bioinformatics_tools` package has the following structure:

 ```python
  class BiologicalSequence(ABC) - allows validate the correctness of the entered sequence,  returns its lenght, symbol by index and right string representation
  class DNASequence(NucleicAcidSequence) and RNASequence(NucleicAcidSequence) - works with nucleotide and amino acid sequences, respectively
  class AminoAcidSequence(BiologicalSequence) - allows to overwrite single letter amino acid sequences to three letter code 
  class FastQFilter - works as filtrator of FastQ-files

```

Example of using `NucleicAcidSequence`, which is acceptable for one or more input sequences:
```
try:
    seq = NucleicAcidSequence("ACtg")
    print(f"Длина последовательности: {len(seq)}")
    print(f"Элемент по индексу 1: {seq[1]}")
    print(f"Строковое представление: {seq}")
    print(seq.complement())
except InvalidSequenceError as e:
    print(e)

###
Длина последовательности: 4
Элемент по индексу 1: C
Строковое представление: ACTG
TGAC
###


try:
    seq = AminoAcidSequence("ACtwtg")
    print(f"Длина последовательности: {len(seq)}")
    print(f"Элемент по индексу 1: {seq[1]}")
    print(f"Строковое представление: {seq}")
    print(seq.translate())
except InvalidSequenceError as e:
    print(e)

### Длина последовательности: 6
Элемент по индексу 1: C
Строковое представление: ACTWTG
AlaCysThrTrpThrGly
###

 The complement of a sequence is given by replacing each base with its complementary base, as given by the dictionary:
    {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'U': 'A', 'u': 'a', 'a': 't', 'c': 'g', 't': 'a', 'g': 'c'}
```

`FastQFilter` module filters a FASTQ file based on GC%, length and quality thresholds and writes the filtered sequences to a new fastq file. It contains such functions as:

- `filter_by_gc_content` takes a list of integers representing GC% borders and filters a fastq file based on those borders. If the list contains one element, it is considered as a high GC% border and 
    all sequences with GC% below it are filtered out. If the list contains two elements, the function filters out all sequences with GC% below the first element and above the second one.

-  `filter_by_length` takes a list of integers representing length borders and filters a fastq file based on those borders. If the list contains one element, it is considered as a low length 
     border and all sequences with length below it are filtered out. If the list contains two elements, the function filters out all sequences with length below the first element and above the second one.

-  `filter_by_quality` takes an integer representing a quality threshold and filters a fastq file based on that threshold. The function calculates the quality of each read in the file, and if 
     the quality is above the given threshold, the read is included in the filtered dictionary.













