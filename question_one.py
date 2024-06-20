from argparse import ArgumentParser
from dataclasses import dataclass, field
from os.path import join, dirname, realpath

import numpy as np


parser = ArgumentParser()
parser.add_argument(
    '-f',
    '--file-path',
    type=str,
    help='The file path to integer sequence input.',
    default=join(dirname(realpath(__file__)), 'data', 'question_one.txt')
)
args = parser.parse_args()


@dataclass
class IntegerSequence:
    '''Class for keeping track of an integer sequence'''
    
    length: int = 0
    sequence: list[int] = field(default_factory=list)
    

class IntegerSequenceReader:
    '''
    Reads an integer sequence file with format:
    >N
    0 1 2 3 4 ... N - 1
    '''
    
    IntegerSequences = list[IntegerSequence]


    class InvalidHeaderError(Exception):
        def __init__(self):
            super().__init__('Header value must be greater than 1.')
            
    class InvalidLengthError(Exception):
        def __init__(self):
            super().__init__('Sequence must have one fewer element than the header value.')


    def __init__(self, filepath: str):
        self.filepath = filepath
        
    def create_integer_sequence(self, header: str) -> IntegerSequence:
        '''Creates and validates an IntegerSequence.'''
        
        # Extract length from header
        len = header.replace('>', '').strip()
        len = int(len)
        
        # Validate length
        if len in [0, 1]:
            raise IntegerSequenceReader.InvalidHeaderError()
        
        return IntegerSequence(length=len)

    def validate_integer_sequences(self, sequences):
        '''Validate all sequence lengths and ensure no element duplication.'''
        
        for seq in sequences:
            if seq.length - 1 != len(seq.sequence):
                raise IntegerSequenceReader.InvalidLengthError()

    def parse(self) -> IntegerSequences:
        '''Parse integer sequence file and initial sequences.'''
        
        try:
            sequences = []
            with open(self.filepath, 'r') as f:
                new_seq = None
                
                for line in f.readlines():
                    if '>' in line:
                        if new_seq:
                            sequences.append(new_seq)
                        
                        new_seq = self.create_integer_sequence(line)
                    else:
                        seq = line.strip()
                        seq = seq.split(' ')
                        new_seq.sequence.extend(map(int, seq))
                        
                sequences.append(new_seq)
                
            self.validate_integer_sequences(sequences)    
            
            return sequences
        except FileNotFoundError:
            print(f'File path: {self.filepath}, does not exist!')


class FindMissingElementsCommand:
    '''Command pattern to run all logic related to question one, also making it testable.'''
    
    def __init__(self, filepath: str):
        self.reader = IntegerSequenceReader(filepath)
    
    def scan_sequence(self, seq: IntegerSequence) -> int:
        '''Scan sequence and return missing value.'''
        
        idx = np.zeros(seq.length, dtype=bool)
        
        for i in seq.sequence:
            idx[i-1] = True
        
        missing_element = np.where(idx == False)[0][0] + 1
            
        return missing_element
    
    def run(self) -> list[int]:
        '''Parse input file and return a list of all missing elements.'''
        
        results = []
        sequences = self.reader.parse()
            
        for int_seq in sequences:
            missing_element = self.scan_sequence(int_seq)
            results.append(missing_element)
            
        return results
            

if __name__ == '__main__':
    command = FindMissingElementsCommand(args.file_path)
    results = command.run()
    
    for result in results:
        print(f'Missing element: {result}')
