import unittest
import os
from filter_logger import FastQFilter


class TestFastQFilter(unittest.TestCase):
    """Test suite for the FastQFilter class."""

    def setUp(self):
        """Create a temporary FASTQ file for testing."""
        self.test_fastq_file = "test.fastq"
        with open(self.test_fastq_file, "w") as f:
            f.write("@SEQ_ID\nACGT\n+\nIIII\n")
            f.write("@SEQ_ID2\nACGTTGCCGA\n+\nIIIIIIII\n")

    def tearDown(self):
        """
        Remove the temporary FASTQ-file created for testing purposes
        after each test to ensure a clean test environment.
        """
        if os.path.exists(self.test_fastq_file):
            os.remove(self.test_fastq_file)

    def test_filter_by_gc_content(self):
        """
        Test the filter_by_gc_content() method of the FastQFilter class.

        The test checks that the method correctly filters sequences based on
        their GC-content. The test uses the temporary FASTQ-file created in
        the setUp() method and checks that the method returns the correct
        number of sequences with GC-content above the specified threshold.
        """
        filter_ = FastQFilter(self.test_fastq_file)

        filtered_records = filter_.filter_by_length(0, 5)
        self.assertEqual(len(filtered_records), 1)

        filtered_records = filter_.filter_by_length(3, 30)
        self.assertEqual(len(filtered_records), 1)

        filtered_records = filter_.filter_by_length(0, 3)
        self.assertEqual(len(filtered_records), 0)


if __name__ == "__main__":
    unittest.main()
