import pytest
from filter_logger import FastQFilter  


@pytest.fixture
def example_fastq_file(tmp_path):
    fastq_content = """@SEQ_ID_1
GATTTGGGGTTTTCCCAGTCACGATC
+
!!!!!!!!!!!!!!!'''''''''''
@SEQ_ID_2
GATTTGGGGTTTTCCCAGTCACG
+
!!!!!!!!!!!!!!!''''''''
@SEQ_ID_3
GATTTGGGGTTTTCCCAGT
+
!!!!!!!!!!!!!!!!!''
"""
    fastq_file = tmp_path / "example.fastq"
    with open(fastq_file, "w") as f:
        f.write(fastq_content)
    return str(fastq_file)


def test_filter_by_length(example_fastq_file):
    """
    Test the filter_by_length() method of the FastQFilter class.

    The test uses a temporary FASTQ file created by the example_fastq_file fixture.
    """

    filter_ = FastQFilter(example_fastq_file)

    filtered_records = filter_.filter_by_length(20, 30)

    assert len(filtered_records) == 2  
    assert filtered_records[0].id == "SEQ_ID_1"  
    assert filtered_records[1].id == "SEQ_ID_2"
    assert filtered_records[1].id != "SEQ_ID_3"
    