from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import GC


class InvalidSequenceError(Exception):
    """Ошибка, если введеная последовательность не соответствует \
        нуклеиновым кислотам"""

    pass


class BiologicalSequence(ABC):
    nucleotides = ["A", "a", "T", "t", "C", "c", "G", "g", "U", "u"]
    amino_acids = [
        "A",
        "R",
        "N",
        "D",
        "C",
        "E",
        "Q",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
    ]

    @abstractmethod
    def __init__(self, sequence: str):
        self.sequence = sequence.upper()
        self.validate_sequence()

    @abstractmethod
    def validate_sequence(self):
        pass

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return f"{self.__class__.__name__}({self.sequence})"


class NucleicAcidSequence(BiologicalSequence):

    def __init__(self, sequence: str):
        super().__init__(sequence)

    def validate_sequence(self):
        for letter in self.sequence:
            if letter not in self.nucleotides:
                raise InvalidSequenceError(
                    f"Некорректный символ {letter} в последовательности"
                )

    def complement(self):
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
        for base in self.sequence:
            new_seq.append(compl_dict[base])
        return "".join(new_seq)

    def reverse(self):
        return self.sequence[::-1]

    def reverse_complement(self):
        seq = self.sequence[::-1]
        return self.complement(seq)


class DNASequence(NucleicAcidSequence):

    def __init__(self, sequence: str):
        super().__init__(sequence)

    def transcribe(self):
        return NucleicAcidSequence.complement(self.sequence)


class RNASequence(NucleicAcidSequence):

    def __init__(self, sequence: str):
        super().__init__(sequence)


class AminoAcidSequence(BiologicalSequence):
    """Реализует перевод однобуквенной АК последовательности в трехбуквенную"""

    one_to_three_letter = {
        "A": "Ala",
        "R": "Arg",
        "N": "Asn",
        "D": "Asp",
        "C": "Cys",
        "E": "Glu",
        "Q": "Gln",
        "G": "Gly",
        "H": "His",
        "I": "Ile",
        "L": "Leu",
        "K": "Lys",
        "M": "Met",
        "F": "Phe",
        "P": "Pro",
        "S": "Ser",
        "T": "Thr",
        "W": "Trp",
        "Y": "Tyr",
        "V": "Val",
    }

    def __init__(self, sequence: str):
        super().__init__(sequence)

    def validate_sequence(self):
        for letter in self.sequence:
            if letter not in self.amino_acids:
                raise InvalidSequenceError(
                    f"Некорректный символ {letter} в последовательности"
                )

    def translate(self):
        three_letter_sequence = []
        for aminoacid in self.sequence:
            three_letter_sequence.append(self.one_to_three_letter[aminoacid])
        return "".join(three_letter_sequence)


class FastQFilter:
    def __init__(self, input_file: str):
        self.input_file = input_file

    def filter_by_length(self, min_length: int, max_length: int):
        filtered_records = []
        for record in SeqIO.parse(self.input_file, "fastq"):
            if min_length <= len(record.seq) <= max_length:
                filtered_records.append(record)
        return filtered_records

    def filter_by_quality(self, min_quality: float):
        filtered_records = []
        for record in SeqIO.parse(self.input_file, "fastq"):
            quality_scores = record.letter_annotations["phred_quality"]
            avg_quality = sum(quality_scores) / len(quality_scores)
            if avg_quality >= min_quality:
                filtered_records.append(record)
        return filtered_records

    def filter_by_gc_content(self, min_gc: float, max_gc: float):
        filtered_records = []
        for record in SeqIO.parse(self.input_file, "fastq"):
            gc_content = GC(record.seq)
            if min_gc <= gc_content <= max_gc:
                filtered_records.append(record)
        return filtered_records

    def save_filtered_records(self, records, output_file: str):
        with open(output_file, "w") as output_handle:
            SeqIO.write(records, output_handle, "fastq")
