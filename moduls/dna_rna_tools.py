def transcribe(*seqs: str) -> str:
    new_seq = []
    sequences = list(seqs)
    for seq in sequences:
        for base in seq:
            rev_seq = base.replace("t", "u").replace("T", "U")
            new_seq.append(rev_seq)
        return new_seq


def reverse(*seqs: str) -> str:
    new_seq = []
    sequences = list(seqs)
    for seq in sequences:
        new = seq[::-1]
        rev_new = " ".join(new)
        new_seq.append(rev_new)
    return new_seq


def complement(*seqs: str) -> str:
    new_seq = []
    compl_dict = {
        "A": "T",
        "C": "G",
        "T": "A",
        "G": "C",
        "U": "A",
        "u": "a",
        "a": "t",
        "c": "g",
        "t": "a",
        "g": "c",
    }
    sequences = list(seqs)
    for seq in sequences:
        for base in seq:
            compl_seq = [compl_dict[base] for base in base]
            compl_seq_new = "".join(compl_seq)
            new_seq.append(compl_seq_new)
        return new_seq


def reverse_complement(*seqs: str) -> str:
    return reverse(complement(seqs))
