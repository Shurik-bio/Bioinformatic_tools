reverse, transcribe, complement, reverse_complement = import(absolute_path = '/home/shurik/HW4/modules').from('reverse', 'transcribe', 'complement', 'reverse_transcribe')
gc_bounds, length_bounds, quality_threshold = import(absolute_path = '/home/shurik/HW4/modules').from('gc_bounds', 'length_bounds', 'quality_threshold')

import string
from os import path
from typing import Union

path_fatsq = "/home/shurik/course_materials/Homeworks/HW5_Files/example_data"
filename = "example_fastq.fastq"
fullparth = os.path.join(path_fatsq, filename)


dir_filt = os.mkdir('filtered')
path_filt = "/home/shurik/course_materials/Homeworks/HW5_Files/example_data"
filename_filt = "filtered_fastq.fastq"
fullpath_to_filt = os.path.join(dir_filt, filename_filt)


def run_dna_rna_tools(*args) -> str:
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

def filter_fastq(input_fastq, output_fastq, length_bounds: Union[int, tuple], gc_bounds: Union[int, tuple], quality_threshold: float) -> dict:
    """
    Filter_fastq is a function which filters a fastq file based on GC%, length and quality thresholds
    and writes the filtered sequences to a new fastq file.

    Parameters
    ----------
    input_fastq : str
        A path to the input fastq file
    output_fastq : str
        A path to the output fastq file
    gc_bounds : int or tuple
        A tuple of integers representing GC% borders
    length_bounds : int or tuple
        A tuple of integers representing length borders
    quality_threshold : float
        A float representing a quality threshold

     length_bounds() function: filters FASTQ records by length
     gc_bounds() function: filters FASTQ records by GC%
     quality_threshold() function: filters FASTQ records by quality

    Returns
    -------
    output_fastq : dict
        A dictionary where the keys are the header lines of the fastq file and
        the values are tuples of sequences and quality strings, filtered based on
        the GC%, length and quality thresholds
    """
    with open(fullparth, 'r', encoding="utf-8") as input_fastq, open(fullpath_to_filt, 'w', encoding="utf-8") as output_fastq:
        fastq_seqs = {}
        name_list = []
        sequence_list = []
        quality_list = []
        for line in input_fastq:
            if line.startswith('@'):
                name_list.append(line)
            elif not line.startswith('+') and not line.startswith('@'):
                sequence_list.append(line)
                quality_list.append(line)
        for key, val1, val2 in zip(name_list, sequence_list, quality_list):
            fastq_seqs[key] = val1, val2
        
        gc_bounds_filtred = gc_bounds(gc_borders)
        length_bounds_filtred = length_bounds(length_borders)
        quality_filtered = quality_threshold(quals)
        
        for key, val1, val2 in zip(name_list, length_bounds_filtred.values(), quality_filtered.values()):
            output_fastq.write('@' + key + '\n' + val1 + '\n' + val2 + '\n')
