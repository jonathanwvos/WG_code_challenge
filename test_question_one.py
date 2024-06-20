from question_one import FindMissingElementsCommand, IntegerSequenceReader
from os.path import join, dirname, realpath

from unittest import TestCase, main


DIR_PATH = dirname(realpath(__file__))
TEST_PATH = join(DIR_PATH, 'test_data')


class TestQuestionOne(TestCase):
    def test_IntegerSequenceReader_ZeroEdgecasesRaisesError(self):
        test_filepath = join(TEST_PATH, 'q1_zero_edgecase.txt')
        command = FindMissingElementsCommand(test_filepath)
        
        with self.assertRaises(IntegerSequenceReader.InvalidHeaderError):
            command.run()
            
    def test_IntegerSequenceReader_DifferentSequenceLengthRaisesError(self):
        test_filepath = join(TEST_PATH, 'q1_invalid_sequence_length.txt')
        command = FindMissingElementsCommand(test_filepath)
        
        with self.assertRaises(IntegerSequenceReader.InvalidLengthError):
            command.run()
            
    def test_IntegerSequenceReader_SucceedsWithValidFile(self):
        test_filepath = join(TEST_PATH, 'q1_valid_file.txt')
        command = FindMissingElementsCommand(test_filepath)
        results = command.run()

        self.assertEqual(results[0], 2)
        self.assertEqual(results[1], 1)
        self.assertEqual(results[2], 146)

        
if __name__ == '__main__':
    main()