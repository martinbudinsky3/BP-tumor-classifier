import pandas as pd


def process_lengths_data(lengths):
    lengths.rename(columns = { 0: 'Chromosome', 1: 'Length' }, inplace = True)
    lengths = lengths.set_index('Chromosome')
    lengths = lengths.loc[:'Y', :'Length']
    
    return lengths

