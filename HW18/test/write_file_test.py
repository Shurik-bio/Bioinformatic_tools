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


def test_save_filtered_records(example_fastq_file, tmp_path):
    """
    Test that save_filtered_records saves records to a file
    """
    
    filter = FastQFilter(example_fastq_file)

    filtered_records = filter.filter_by_length(26, 26)
    output_file = tmp_path / "filtered_output.fastq"

    filter.save_filtered_records(filtered_records, str(output_file))

    assert output_file.is_file()

    records = list(SeqIO.parse(str(output_file), "fastq"))
    assert len(records) == 2
    assert records[0].id == "SEQ_ID_1"
    assert str(records[0].seq) == "GATTTGGGGTTTTCCCAGTCACGATC"
    assert records[1].id == "SEQ_ID_2"
    assert str(records[1].seq) == "GATTTGGGGTTTTCCCAGTCACGATC"
