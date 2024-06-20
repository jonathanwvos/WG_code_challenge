from question_one import IntegerSequenceReader
from os.path import join, dirname, realpath

from unittest import TestCase, main


DIR_PATH = dirname(realpath(__file__))
TEST_PATH = join(DIR_PATH, 'test_data')


class TestQuestionOne(TestCase):
    def test_IntegerSequenceReader_ZeroEdgecasesRaisesError(self):
        
        with self.assertRaises(IntegerSequenceReader.InvalidHeaderError):
            test_filepath = join(TEST_PATH, 'q1_zero_edgecase.txt')
            IntegerSequenceReader(test_filepath)
            
    def test_IntegerSequenceReader_DifferentSequenceLengthRaisesError(self):
        with self.assertRaises(IntegerSequenceReader.InvalidLengthError):
            test_filepath = join(TEST_PATH, 'q1_invalid_sequence_length.txt')
            IntegerSequenceReader(test_filepath)
            
    def test_IntegerSequenceReader_SucceedsWithValidFile(self):
        raise NotImplementedError

        
if __name__ == '__main__':
    main()