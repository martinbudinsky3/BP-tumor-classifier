import pandas as pd

from gap_data_processor import GapDataProcessor
from lengths_data_processor import process_lengths_data

Mb = 1000000
LENGTHS_PATH = '../datasets/hs37d5.fa.fai'
GAP_PATH = '../datasets/gap.txt'

def init_lengths():
    lengths = pd.read_csv(LENGTHS_PATH, sep='\t', header=None)
    return process_lengths_data(lengths)
    
    
def init_centromeres():
    gdp = GapDataProcessor(GAP_PATH)
    return gdp.get_centromeres()


def init_chromosome_names():
    chromosome_names = [str(_chr) for _chr in range(1, 23)]
    chromosome_names.append('X')

    return chromosome_names