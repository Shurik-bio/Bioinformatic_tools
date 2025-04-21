import pytest
from Bio import SeqIO
from filter_logger import FastQFilter  # Замените на имя вашего модуля


@pytest.fixture
def example_fastq_file(tmp_path):
    # Создаем временный FASTQ файл для тестирования
    """
    Fixture, returning a path to a temporary FASTQ-file with two sequences
    The file is deleted after the test.
    """
    fastq_content = """@SEQ_ID_1
GATTTGGGGTTTTCCCAGTCACGATC
+
!''*((((***+))%%%++)(*****
@SEQ_ID_2
GATTTGGGGTTTTCCCAGTCACGATC
+
!''*((((***+))%%%++)(*****
"""
    fastq_file = tmp_path / "example.fastq"
    with open(fastq_file, "w") as f:
        f.write(fastq_content)
    return str(fastq_file)


def test_read_fastq(example_fastq_file):
    """
    Test that read_fastq reads sequences from a file
    """
    
    filter = FastQFilter(example_fastq_file)

    records = list(SeqIO.parse(filter.input_file, "fastq"))

    assert len(records) == 2
    assert records[0].id == "SEQ_ID_1"
    assert str(records[0].seq) == "GATTTGGGGTTTTCCCAGTCACGATC"
    assert records[1].id == "SEQ_ID_2"
    assert str(records[1].seq) == "GATTTGGGGTTTTCCCAGTCACGATC"
