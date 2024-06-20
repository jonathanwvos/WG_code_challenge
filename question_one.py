from argparse import ArgumentParser
from dataclasses import dataclass, field

parser = ArgumentParser()
parser.add_argument(
    '-f',
    '--file-path',
    type=str,
    help='The file path to integer sequence input.'
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
    
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.sequences = []
        self.parse_file()

    def parse_file(self):
        try:
            with open(self.filepath, 'r') as f:
                new_seq = None
                
                for line in f.readlines():
                    if '>' in line:
                        if new_seq:
                            self.sequences.append(new_seq)
                        
                        new_seq = IntegerSequence()
                        
                        len = line.replace('>', '').strip()
                        new_seq.length = int(len)
                    else:
                        seq = line.strip()
                        seq = seq.split(' ')
                        new_seq.sequence.extend(map(int, seq))
                        
                self.sequences.append(new_seq)
                    
        except FileNotFoundError:
            print(f'File path: {self.filepath}, does not exist!')


if __name__ == '__main__':
    reader = IntegerSequenceReader(args.file_path)
