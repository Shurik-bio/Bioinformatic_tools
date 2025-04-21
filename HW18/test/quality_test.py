import unittest
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from filter_logger import FastQFilter


class TestFastQFilter(unittest.TestCase):

    def setUp(self):
        # Создаем временный файл для тестирования
        """
        Create a temporary file for testing. The file contains two sequences in FASTQ format
        with different quality scores. The file is deleted after the test.
        """

        self.test_fastq_file = "test.fastq"
        with open(self.test_fastq_file, "w") as f:
            f.write(
                """@SEQ_ID_1
GATTTGGGGTTTAAAGGGTTT
+SEQ_ID_1
!''*((((***+))%%%++((
@SEQ_ID_2
GATTTGGGGTTTAAAGGGTTT
+SEQ_ID_2
!''*((((***+))%%%
"""
            )

    def tearDown(self):
        """
        Remove the temporary FASTQ file created for testing purposes
        after each test to ensure a clean test environment.
        """

        if os.path.exists(self.test_fastq_file):
            os.remove(self.test_fastq_file)

    def test_filter_by_gc_content(self):
        """
        Test the filter_by_quality() method of the FastQFilter class.

        The test checks that the method correctly filters sequences based on
        their quality scores. The test uses the temporary FASTQ file created
        in the setUp() method and checks that the method returns the correct
        number of sequences with quality scores above the specified threshold.
        """

        filter = FastQFilter(self.test_fastq_file)

        filtered_records = filter.filter_by_quality(min_quality=3)
        self.assertEqual(len(filtered_records), 1)

        filtered_records = filter.filter_by_quality(min_quality=7)
        self.assertEqual(len(filtered_records), 1)

        filtered_records = filter.filter_by_quality(min_quality=40)
        self.assertEqual(len(filtered_records), 0)


if __name__ == "__main__":
    unittest.main()
