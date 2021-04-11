from utils import init_lengths, init_chromosome_names, Mb

lengths = init_lengths()
chromosome_names = init_chromosome_names()


def loh(data, LOH_TRESHOLD=15*Mb):
    long_lohs = 0
    for _chr in chromosome_names:
        chr_data = data.loc[data['Chromosome'] == _chr]
        long_loh_segments = chr_data.loc[(chr_data['Length'] > LOH_TRESHOLD) & (chr_data['Length'] < lengths.loc[_chr, 'Length'])]
        long_lohs += len(long_loh_segments.index)

    return long_lohs