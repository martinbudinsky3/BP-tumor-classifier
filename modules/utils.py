import pandas as pd
from pathlib import Path

from .gap_data_processor import GapDataProcessor
from .lengths_data_processor import process_lengths_data

here = Path(__file__).parent
LENGTHS_PATH = here/'data/hs37d5.fa.fai'
GAP_PATH = here/'data/gap.txt'
Mb = 1000000


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