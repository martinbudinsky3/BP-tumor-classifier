from .utils import init_lengths, init_chromosome_names, Mb

lengths = init_lengths()
chromosome_names = init_chromosome_names()


def loh(data, LOH_TRESHOLD=15*Mb):
    """
    Implementation of LOH method
    
    Parameters
    ----------
    data: pandas.DataFrame
        DataFrame containing preprocessed segmental report data
        
    LOH_TRESHOLD=15000000: int, optional
        Value of parameter LOH_TRESHOLD (treshold length for segments that count to LOH score) of LOH method.
        If not provided, LOH is count for default value of parameter 15 Mb
    
    Returns
    -------
    loh_score: int
        LOH score of sample
    """
    
    long_lohs = 0
    for _chr in chromosome_names:
        chr_data = data.loc[data['Chromosome'] == _chr]
        long_loh_segments = chr_data.loc[(chr_data['Length'] > LOH_TRESHOLD) & (chr_data['Length'] < lengths.loc[_chr, 'Length'])]
        long_lohs += len(long_loh_segments.index)

    return long_lohs