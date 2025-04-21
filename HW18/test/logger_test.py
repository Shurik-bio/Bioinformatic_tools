import pytest
import logging
from filter_logger import FastQFilter 

@pytest.fixture
def example_fastq_file(tmp_path):
    """
    Fixture that creates a temporary FASTQ-file for testing.

    The file has two sequences with the same identifier and the same sequence
    but different quality scores.

    Returns the path to the temporary file as a string.
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

def test_logging(caplog, example_fastq_file):
    """
    Test that the filter_by_length method logs a message.

    The method is called with a valid FASTQ file and a length range.
    The caplog fixture is used to capture the log messages.
    The length of the log messages is asserted to be greater than 0.
    """
    
    filter = FastQFilter(example_fastq_file)

    with caplog.at_level(logging.INFO):
        filter.filter_by_length(10, 30) 

    assert len(caplog.records) > 0 
    
    
