from typing import Union

def gc_bounds(gc_borders: Union[int, tuple]) -> dict:
    if len(gc_borders) == 1:
        gc_borders = (0, gc_borders[0])
        gc_bounds_filtred = {}
        for read in fastq_seqs.values():
            reads = str(read[0])
            gc_content = float(
                ((reads.count("C") + reads.count("G")) / len(reads)) * 100
            )
            if gc_borders[0] <= gc_content <= gc_borders[1]:
                gc_bounds_filtred = {
                    key: val for key, val in fastq_seqs.items() if val == read
                }
        return gc_bounds_filtred


def length_bounds(length_borders: Union[int, tuple]) -> dict:
    if len(length_borders) == 1:
        length_borders = (0, length_borders[0])
        length_bounds_filtred = {}
        for read in gc_bounds_filtred.values():
            reads = str(read[0])
            length = len(reads)
            if length_borders[0] <= length <= length_borders[1]:
                length_gc_bounds_filtred = {
                    key: val for key, val in gc_bounds_filtred.items() if val == read
                }
        return length_gc_bounds_filtred


def quality_threshold(quals: float) -> dict:
    for read_quals in length_gc_bounds_filtred.values():
        quality_counter = 0
        reads = str(read_quals[1])
        for symbol in reads:
            quality = ord(symbol) - 33
            quality_counter += quality
        quality_read = quality_counter / len(reads)
        if quality_read >= quals:
            quality_filtered = {
                key: val
                for key, val in length_gc_bounds_filtred.items()
                if val[1] == reads
            }
    return quality_filtered
