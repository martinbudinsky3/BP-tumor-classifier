from .utils import init_lengths, init_centromeres, init_chromosome_names, Mb

lengths = init_lengths()
centromeres = init_centromeres()
chromosome_names = init_chromosome_names()


def loh(data, LOH_TRESHOLD=15*Mb, with_centromere=True):
    """
    Implementation of LOH method
    
    Parameters
    ----------
    data: pandas.DataFrame
        DataFrame containing preprocessed segmental report data
        
    LOH_TRESHOLD=15000000: int, optional
        Value of parameter LOH_TRESHOLD (treshold length for segments that count to LOH score) of LOH method.
        If not provided, LOH is count for default value of parameter 15 Mb
        
    with_centromere=True: bool, optional
        Flag indicating the way in which condition, that only segments shorter than chromosome length are counted to result score, should be applied.
        If the flag is set to True, algorithm will not count to result score only segments with length of whole chromosome (even with centromere),
        assuming that variant caller joins segments adjacent with centromere.
        Although if the flag is set to False, algorithm will not count to result score even segments their lenghts sums to length of chromosome
        (without centromere).
        
    
    Returns
    -------
    loh_score: int
        LOH score of sample
    """
    
    long_lohs = 0
    for _chr in chromosome_names:
        chr_data = data.loc[data['Chromosome'] == _chr]
        chr_len = lengths.loc[_chr, 'Length']
        
        long_loh_segments = chr_data.loc[(chr_data['Length'] > LOH_TRESHOLD)]
        if with_centromere:
            chr_long_loh_count = count_chr_long_lohs_with_centromere(long_loh_segments, chr_len)
        else:
            chr_long_loh_count = count_chr_long_lohs_without_centromere(long_loh_segments, _chr, chr_len, centromeres)
                
        long_lohs += chr_long_loh_count
            
        
    return long_lohs


def count_chr_long_lohs_with_centromere(long_loh_segments, chr_len):
    long_loh_segments = long_loh_segments.loc[long_loh_segments['Length'] < chr_len]

    return len(long_loh_segments.index)
    
    
def count_chr_long_lohs_without_centromere(long_loh_segments, _chr, chr_len, centromeres):
    centromere_start = centromeres.loc[centromeres['Chromosome'] == _chr, 'Start'].iloc[0]
    centromere_end = centromeres.loc[centromeres['Chromosome'] == _chr, 'End'].iloc[0]
    centromere_len = centromere_end - centromere_start
    chr_len_without_centromere = chr_len - centromere_len

    sum_of_long_loh_segments_lengths = long_loh_segments['Length'].sum()

    if sum_of_long_loh_segments_lengths < chr_len_without_centromere:
        return len(long_loh_segments.index)
    
    return 0