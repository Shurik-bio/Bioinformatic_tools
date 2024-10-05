def gc_bounds(x, y: float) -> dict:
    for read in seqs.values():
        reads = str(read[0])
    	gc_content = float(((reads.count("C") + reads.count("G"))/len(reads))*100)
        if gc_content <= x or x <= gc_content <= y:
            gc_bounds_filtred = {key: val for key, val in seqs.items() if val == read}
        elif x == 0:
            gc_bounds_filtred = {key: val for key, val in seqs.items() if val == read}
        return gc_bounds_filtred

def length_bounds(x, y: int) -> dict:
    for read in seqs.values():
        reads = str(read[0])
        len_reads = len(reads)
        if len_reads <= x or x <= len_reads <= y:
            length_filtered = {key: val for key, val in seqs.items() if val == read}
        elif x == 0:
            length_filtered = {key: val for key, val in seqs.items() if val == read}
        return length_filtered


def quality_threshold(x: float) -> dict:
    for read_quals in seqs.values():
        quality_counter = 0
        reads = str(read_quals[1])
        for symbol in reads:
            qual = ord(symbol) - 33
            quality_counter += qual
            quality_read = quality_counter/len(reads)
        if quality_read >= x:
            quality_filtered = {key: val for key, val in seqs.items() if val[1] == reads}
        return quality_filtered
