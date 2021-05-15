import pandas as pd


def process_lengths_data(lengths):
    """
    Function for preprocessing of chromosomes lengths data
    
    Parameters
    ----------
    lenths : pandas.DataFrame
        DataFrame with unprocessed lengths data
        
    Returns
    -------
    lengths: pandas.DataFrame
        DataFrame with preprocessed lengths data
    """
    
    lengths.rename(columns = { 0: 'Chromosome', 1: 'Length' }, inplace = True)
    lengths = lengths.set_index('Chromosome')
    lengths = lengths.loc[:'Y', :'Length']
    
    return lengths

