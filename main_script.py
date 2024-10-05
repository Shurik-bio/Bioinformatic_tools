from moduls.dna_rna_tools import reverse, transcribe, complement, reverse_complement
from moduls.filter_fastq import gc_bounds, length_bounds, quality_threshold
import string

def run_dna_rna_tools(*args):
"""
    Applies a given function from the dna_rna_tools module to a variable number of sequences.

    Args:
        seqs (str): One or more DNA/RNA sequences to be processed
        function (str): The function to be applied to the sequences. Must be one of 'transcribe', 'complement', 'reverse', or 'reverse_complement'

    Returns:
        str: The processed sequence(s). If only one sequence is provided, it is returned as a single string. If multiple sequences are provided, they are returned as a string with the sequences separated by spaces.

    The following functions are currently supported:
        - Transcribe DNA sequences to RNA.

    Given one or more DNA sequences, return the equivalent RNA sequences.

  Parameters:
  *seqs : str
    one or more DNA sequences

  Returns:
  a list of RNA sequences

  Examples:
  >>> transcribe('ATGC')
  ['AUGC']
  >>> transcribe('ATGC', 'ATGC')
  ['AUGC', 'AUGC']

        - Reverse a DNA or RNA sequences.

    Parameters:
    seqs (str): The DNA/RNA sequence(s) to reverse.

    Returns:
    str: The reversed DNA sequence(s), as a string.

    Example:
    >>> reverse('ATCG')
    'GCAT'
    >>> reverse('ATCG', 'TGCA')
    ['GCAT', 'ACTG']

        - Returns the complement of one or more sequences.

    The complement of a sequence is given by replacing each base with its
    complementary base, as given by the dictionary

    {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'U': 'A', 'u': 'a', 'a': 't', 'c': 'g', 't': 'a', 'g': 'c'}

    Parameters:
    one or more DNA/RNA sequences

    Returns:
    list of complement sequences

        - Returns the reverse complement of one or more sequences.

    The reverse complement of a sequence is given by reversing the sequence
    and then taking the complement of the reversed sequence.

    Parameters:
    one or more DNA/RNA sequences

    Returns:
    list of reverse complement sequences

  """
    *seqs, function = args
        if function == 'transcribe':
            result = transcribe(seqs)
        elif function == 'complement':
            result = complement(seqs)
        elif function == 'reverse_complement':
            result = reverse_complement(seqs)
        elif function == 'reverse':
            result = reverse(seqs)
        if len(result) == 1:
            return ''.join(result)
        else:
            return ' '.join(result)

def filter_fastq(seqs: dict, length_bounds: int, gc_bounds: float, quality_threshold: float) -> dict:
    """
    Filter FASTQ records by one or more criteria.

    Each criterion is given as a function, which takes a single FASTQ record
    as argument and returns True or False. The functions are applied
    sequentially, and a record is only passed through if all criteria return
    True.

    The following criteria are currently supported:

    - `length_bounds`: a tuple of two integers, representing the minimum and
      maximum length of the read

      length_bounds() function: filters FASTQ records by length

  Args:
      x (int): the minimum length of the read
      y (int): the maximum length of the read

  Returns:
      dict: Filtered FASTQ records by length

    - `seqs`: a list of strings, representing sequence motifs to be searched
      for
    - `gc_bounds`: a tuple of two floats, representing the minimum and maximum

      gc_bounds() function: calculates GC-content of the read and filters according to a given interval

  Args:
      x (float): the minimum GC content
      y (float): the maximum GC content

  Returns:
      dict: Filtered FASTQ records 

    - `quality_threshold`: a float, representing the minimum average basecall
      quality of the read

      Filters FASTQ records by quality average basecall quality of a read.
    Receive the integer quality scores of a read and return the average quality for that read.

  Args:
      x (float): the minimum average basecall quality

  Returns:
      dict: filtered FASTQ records with given quality threshold

    Returns a generator of filtered FASTQ records.
    """
  passed_reads : dict = {}
  for read in seqs.values():
    sequence = length_bounds, gc_bounds
    quality = quality_threshold
    if sequence == length_bounds(sequence):
      return length_bounds_filtred
    elif sequence == gc_bounds(sequence):
      return gc_bounds_filtred
    elif quality == quality_threshold(quality):
      return quality_filtered
  for key, value in gc_bounds_filtred.items():
    for key, value in length_filtered.items():
      for key, value in quality_filtered.items():
        if gc_bounds_filtred[key] == length_filtered[key] == quality_filtered[key]:
          passed_reads[key] = value
      return passed_reads
