from .utils import init_lengths, init_centromeres, init_chromosome_names, Mb

lengths = init_lengths()
centromeres = init_centromeres()
chromosome_names = init_chromosome_names()

def tai(data, TELOMERE_SIZE=2*Mb):
    ntai = 0
    for _chr in chromosome_names:
        chr_data = data.loc[data['Chromosome'] == _chr]
        chr_len = lengths.loc[_chr, 'Length']
        centromere_start = centromeres.loc[centromeres['Chromosome'] == _chr, 'Start'].iloc[0]
        centromere_end = centromeres.loc[centromeres['Chromosome'] == _chr, 'End'].iloc[0]
        
        tais_segments = chr_data.loc[
            (chr_data['Start'] < TELOMERE_SIZE) & (chr_data['End'] <= centromere_start) | 
            (chr_data['End'] > chr_len - TELOMERE_SIZE) & (chr_data['Start'] >= centromere_end)
        ]
        ntai += len(tais_segments.index)

    return ntai