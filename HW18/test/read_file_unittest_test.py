import unittest
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from filter_logger import FastQFilter


class TestFastQFilter(unittest.TestCase):

    def setUp(self):
        """
        Create a temporary FASTQ file for testing.
        The file is deleted after the test.
        """
        self.test_fastq_file = "test.fastq"
        with open(self.test_fastq_file, "w") as f:
            f.write("@SEQ_ID\nACGT\n+\nIIII\n")

    def tearDown(self):
        """
        Remove the temporary FASTQ-file created for testing purposes
        after each test to ensure a clean test environment.
        """

        if os.path.exists(self.test_fastq_file):
            os.remove(self.test_fastq_file)

    def test_filter_by_length(self):
        """
        Test the filter_by_length() method of the FastQFilter class.
        """
        filter = FastQFilter(self.test_fastq_file)

        filtered_records = filter.filter_by_length(1, 4)
        self.assertEqual(len(filtered_records), 1)

        filtered_records = filter.filter_by_length(5, 10)
        self.assertEqual(len(filtered_records), 0)


if __name__ == "__main__":
    unittest.main()
