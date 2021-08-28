from .utils import init_lengths, init_centromeres, init_chromosome_names, Mb

lengths = init_lengths()
centromeres = init_centromeres()
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
        chr_len = lengths.loc[_chr, 'Length']
        centromere_start = centromeres.loc[centromeres['Chromosome'] == _chr, 'Start'].iloc[0]
        centromere_end = centromeres.loc[centromeres['Chromosome'] == _chr, 'End'].iloc[0]
        centromere_len = centromere_end - centromere_start
        chr_len_without_centromere = chr_len - centromere_len
        
        long_loh_segments = chr_data.loc[(chr_data['Length'] > LOH_TRESHOLD)]
        sum_of_long_loh_segments_lengths = long_loh_segments['Length'].sum()

        if sum_of_long_loh_segments_lengths < chr_len_without_centromere:
            long_lohs += len(long_loh_segments.index)

    return long_lohs