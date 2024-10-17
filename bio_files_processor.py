import os

path_fasta = "/home/shurik/course_materials/Homeworks/HW5_Files/example_data"
file_fasta_name = "example_multiline_fasta.fasta"
fullparth_fasta = os.path.join(path_fasta, file_name)

path_blast = "/home/shurik/course_materials/Homeworks/HW5_Files/example_data"
file_blast_name = "example_blast_results.txt"
fullparth_blast = os.path.join(path_blast, file_blast_name)


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str) -> str:
    """
    Converts a multi-line fasta file to a one-line fasta file.

    Parameters
    ----------
    input_fasta : str
        A path to the input fasta file
    output_fasta : str
        A path to the output fasta file
    """

    with open(fullparth, "r") as input_fasta, open("output.fasta", "w") as output_fasta:
        sequence = ""
        for line in input_fasta:
            line = line.strip()
            if line.startswith(">"):
                output_file.write(line + "\n")
                sequence = ""
            else:
                sequence += line
        output_fasta.write(sequence)


def parse_blast_output(input_file: str, output_file: str) -> str:
    """
    Parses a BLAST output file and writes the results to an output file.

    Parameters
    ----------
    input_file : str
        A path to the input BLAST output file
    output_file : str
        A path to the output file
    """

    with open(fullparth_blast, "r") as input_file, open(
        "output.txt", "w"
    ) as output_file:
        lines = input_file.readlines()
        protein_describtion = []

        for line in lines:
            if line.startswith("Sequences producing significant alignments:"):
                protein_describtion.append(line.strip())

            if line.strip() == "":
                break

    protein_describtion = sorted(set(protein_describtion))

    output_file.write(protein_describtion)
