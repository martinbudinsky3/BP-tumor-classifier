import pandas as pd
import numpy as np
import scipy.stats as stats

from .utils import init_lengths, init_centromeres, init_chromosome_names, Mb

LST_SMbs = [x for x in range(3, 12)]
lengths = init_lengths()
centromeres = init_centromeres()
chromosome_names = init_chromosome_names()


def lst(data, vcf_reader=None, sample_name=None, LST_SMb_param=11):
    """
    Implementation of LST method
    
    Parameters
    ----------
    data: pandas.DataFrame
        DataFrame containing preprocessed segmental report data
        
    vcf_reader: vcf.Reader, optional
        Instance of vcf.Reader for VCF file of input sample. If not provided, LST is count only based on copy numbers.
        
    sample_name: str, optional
        Name of sample in VCF file. If not provided, LST is count only based on copy numbers.
        
    LST_SMb_param=11: int, optional
        Value of parameter LST_SMb (in Mb) of LST method. If not provided, LST is count for value of parameter 3 - 11 Mb
    
    Returns
    -------
    lst, dna_index: tuple
        If value of LST_SMb_param is provided function returns LST score for provided parameter value
        and DNA index of sample
        
    lsts, dna_index: tuple
        If value of LST_SMb_param is not provided function returns dictionary with LST scores for parameter value 3 - 11 Mb 
        and DNA index of sample. Dictionary with LST scores contains keys in format: LST_<LST_SMb>Mb.
    """
    
    data = fill_segments(data)
    dna_index = count_dna_index(data)
        
    # if sample has vcf data coerce by copy numbers and allelic frequencies else only by copy numbers
    if not vcf_reader is None:
        data = count_allelic_freqs(data, vcf_reader, sample_name)
        data = coercing(data)
    else:
        data = coercing(data, count_af=False)
    
    # count lst only for LST_SMb_param size
    if not LST_SMb_param is None:
        lst = count_lsts(data, LST_SMb_param*Mb)
        return lst, dna_index
    
    # count lst for sizes 3,4...11Mb
    else:
        lsts = {}
        for LST_SMb in LST_SMbs:
            lsts['LST_' + str(LST_SMb)+'Mb'] = count_lsts(data, LST_SMb*Mb)

        return lsts, dna_index      

    
# fill gaps in segmented genome profile with segments with copy number 2
def fill_segments(data):
    filled_data = data.copy()
    index_filled_data = 0
    normal_cn = 2
    
    for index, segment in data.iterrows():
        # first cnv region in chromosome
        if index == 0 or data.loc[index-1, 'Chromosome'] != data.loc[index, 'Chromosome']:
            if segment['Start'] != 0:
                filled_data = insert_row(filled_data, segment['Chromosome'], normal_cn, 0, segment['Start'], index_filled_data)
                index_filled_data += 1

        # not first cnv region in chromosome 
        elif data.loc[index-1, 'End'] != data.loc[index, 'Start']:
            prev = data.loc[index-1]
            filled_data = insert_row(filled_data, segment['Chromosome'], normal_cn, prev['End'], segment['Start'], index_filled_data)
            index_filled_data += 1

        # last cnv region in chromosome
        if index == len(data) - 1 or data.loc[index+1, 'Chromosome'] != data.loc[index, 'Chromosome']:
            chr_len = lengths.loc[segment['Chromosome'], 'Length']
            if segment['End'] != chr_len:      
                filled_data = insert_row(filled_data, segment['Chromosome'], normal_cn, segment['End'], chr_len, index_filled_data + 1)
                index_filled_data += 1

        index_filled_data += 1

    filled_data = remove_centromeres(filled_data)

    return filled_data


# insert new segment into dataframe
def insert_row(data, _chr, cn, start, end, index):
    new_segment = pd.DataFrame({
        'Chromosome': [ _chr ],
        'Copy Number': [ cn ],
        'Start': [ start ],
        'End': [ end ]
    })

    return pd.concat([data.iloc[:index], new_segment, data.iloc[index:]]).reset_index(drop=True)


# remove centromeric regions of chromosomes from segmented profile
def remove_centromeres(data):
    for _chr in chromosome_names:
        centromere_start = centromeres.loc[centromeres['Chromosome'] == _chr, 'Start'].iloc[0]
        centromere_end = centromeres.loc[centromeres['Chromosome'] == _chr, 'End'].iloc[0]

        data = drop_segments_inside_centromere(data, _chr, centromere_start, centromere_end)
        data = cut_segments_overlapping_centromere(data, _chr, centromere_start, centromere_end)
        data = cut_segments_ending_in_centromere(data, _chr, centromere_start, centromere_end)
        data = cut_segments_starting_in_centromere(data, _chr, centromere_start, centromere_end)
        
    data = name_chr_arms(data)

    return data


def drop_segments_inside_centromere(data, _chr, centromere_start, centromere_end):
    inside_centromere_cond = (data['Chromosome'] == _chr) & (data['Start'] >= centromere_start) & (data['End'] <= centromere_end)
    data = data.drop(data[inside_centromere_cond].index).reset_index(drop=True)
    
    return data


def cut_segments_overlapping_centromere(data, _chr, centromere_start, centromere_end):
    overlaps_centromere_from_both_sides_cond = (data['Chromosome'] == _chr) & (data['Start'] < centromere_start) & (data['End'] > centromere_end)
    segments_overlaping_centromere = data.loc[overlaps_centromere_from_both_sides_cond]
    
    if not segments_overlaping_centromere.empty:
        segment_overlaping_centromere = segments_overlaping_centromere.iloc[0]
        data.loc[segment_overlaping_centromere.name, 'End'] = centromere_start
        data = insert_row(
            data, 
            _chr, 
            segment_overlaping_centromere['Copy Number'], 
            centromere_end, 
            segment_overlaping_centromere['End'], 
            segment_overlaping_centromere.name + 1
        )
    
    return data


def cut_segments_ending_in_centromere(data, _chr, centromere_start, centromere_end):
    end_in_centromere_cond = (data['Chromosome'] == _chr) & (data['End'] > centromere_start) & (data['End'] <= centromere_end)
    data.loc[end_in_centromere_cond, 'End'] = centromere_start
    
    return data


def cut_segments_starting_in_centromere(data, _chr, centromere_start, centromere_end):
    start_in_centromere_cond = (data['Chromosome'] == _chr) & (data['Start'] >= centromere_start) & (data['Start'] < centromere_end)
    data.loc[start_in_centromere_cond, 'Start'] = centromere_end
    
    return data


# name chromosome arms with 'p' and 'q' label
def name_chr_arms(data):
    data['Arm'] = 'p'
    for _chr in chromosome_names:
        centromere_end = centromeres.loc[centromeres['Chromosome'] == _chr, 'End'].iloc[0]

        data.loc[(data['Chromosome'] == _chr) & (data['Start'] >= centromere_end), 'Arm'] = 'q'

    return data


# count metric DNA index for sample as average_copy_number / 2
def count_dna_index(data):
    data = update_segments_lengths(data)
    cns = list(data['Copy Number'])
    weights = list(data['Length'])
    normal_cn = 2
    
    for _chr in chromosome_names:
        chr_segments = data[data['Chromosome'] == _chr]
        if chr_segments.empty:
            centromere_start = centromeres.loc[centromeres['Chromosome'] == _chr, 'Start'].iloc[0]
            centromere_end = centromeres.loc[centromeres['Chromosome'] == _chr, 'End'].iloc[0]
            chromosome_end = lengths.loc[_chr, 'Length']

            cns.extend([normal_cn, normal_cn])
            weights.extend([centromere_start, chromosome_end - centromere_end])

    cn_avg = np.average(cns, weights=weights)
    dna_index = cn_avg / 2

    return dna_index


def update_segments_lengths(data):
    data['Length'] = data['End'] - data['Start']
    
    return data


# function for counting allelic frequencies for each segment
def count_allelic_freqs(data, vcf_reader, sample, qual_threshold = 50):
    data_with_af = data.copy()
    data_with_af['Allelic Frequencies'] = [list() for x in range(len(data_with_af.index))]

    for index, segment in data_with_af.iterrows():
        try:
            segment_records = vcf_reader.fetch(segment['Chromosome'], segment['Start'], segment['End'])
            allelic_freqs = []
        
            for record in segment_records:
                sample_data = record.genotype(sample).data

                if has_quality(record) and sample_data.GT != './.' and sample_data.GT != '0/0' and sample_data.AD != './.' and sample_data.AD != None \
                    and sample_data.AD[1] != 0:

                    allelic_freq = sample_data.AD[1] / (sample_data.AD[0] + sample_data.AD[1])
                    allelic_freqs.append(allelic_freq)

            data_with_af.at[index, 'Allelic Frequencies'] = allelic_freqs

        except (ValueError, AttributeError) as e:
            continue

    return data_with_af


# function that checks if vcf record has required quality
def has_quality(record, qual_threshold = 50):
    info = record.INFO
    return record.QUAL > qual_threshold and ('QD' not in info.keys() or info['QD'] > 10.0) and ('MQ' not in info.keys() or info['MQ'] > 40.0) \
        and ('FS' not in info.keys() or info['FS'] < 30.0 ) and ('SOR' not in info.keys() or info['SOR'] < 3.0) \
        and ('MQRankSum' not in info.keys() or info['MQRankSum'] > -12.5) and ('ReadPosRankSum' not in info.keys() or info['ReadPosRankSum'] > -8.0)


# coercing function
def coercing(data, count_af=True, S_small=3*Mb):
    data = update_segments_lengths(data)
    
    while len(data) > 0:
        # get smallest segment
        smallest_segment = data[data['Length'] == data['Length'].min()].iloc[0]
        index = smallest_segment.name

        # filter out?
        if smallest_segment['Length'] < S_small:
            # not first or last segment of profile?
            if index != 0 and index != len(data) - 1:
                prev = data.loc[ index-1 ]
                _next = data.loc[ index+1 ]

                # can link?
                if prev['Chromosome'] == _next['Chromosome'] and prev['Arm'] == _next['Arm'] and prev['Copy Number'] == _next['Copy Number']:
                    # if sample has vcf data check allelic frequencies else join only based on copy number
                    if not count_af:
                        data = link_segments_without_af(data, prev, _next)
                    elif have_equal_allelic_freqs(prev, _next):
                        data = link_segments_with_af(data, prev, _next, smallest_segment)

            # delete small segment
            data = data.drop(index=index).reset_index(drop=True)

        # if there are no small segments left -> end
        else:
            break
        
    return data


# linking adjacent segments of small filtered out segment - segments without allelic frequencies 
def link_segments_without_af(data, prev, _next):
    data = data.drop(index=prev.name)
    data = data.drop(index=_next.name)

    data = insert_row_without_af(data, prev['Chromosome'], prev['Copy Number'],  prev['Start'], _next['End'], \
                            prev['Arm'], prev.name)

    return data


# insert new segment that was created by linking - segment without allelic frequencies
def insert_row_without_af(data, _chr, cn, start, end, arm, index):
    new_segment = pd.DataFrame({
        'Chromosome': [ _chr ],
        'Copy Number': [ cn ],
        'Length': [ end - start ],
        'Start': [ start ],
        'End': [ end ],
        'Arm': [arm]
    })

    return pd.concat([data.iloc[:index], new_segment, data.iloc[index:]]).reset_index(drop=True)


# linking adjacent segments of small filtered out segment 
def link_segments_with_af(data, prev, _next, small):
    data = data.drop(index=prev.name)
    data = data.drop(index=_next.name)

    new_allelic_freqs = []
    new_allelic_freqs.extend(prev['Allelic Frequencies'])
    new_allelic_freqs.extend(_next['Allelic Frequencies'])
    new_allelic_freqs.extend(small['Allelic Frequencies'])

    data = insert_row_with_af(data, prev['Chromosome'], prev['Copy Number'],  prev['Start'], _next['End'], new_allelic_freqs, \
                            prev['Arm'], prev.name)

    return data


# insert new segment that was created by linking 
def insert_row_with_af(data, _chr, cn, start, end, allelic_freqs, arm, index):
    new_segment = pd.DataFrame({
        'Chromosome': [ _chr ],
        'Copy Number': [ cn ],
        'Length': [ end - start ],
        'Start': [ start ],
        'End': [ end ],
        'Allelic Frequencies': [allelic_freqs],
        'Arm': [arm]
    })

    return pd.concat([data.iloc[:index], new_segment, data.iloc[index:]]).reset_index(drop=True)


# statistical tests for equality of allelic frequencies of two segments
def have_equal_allelic_freqs(segment1, segment2):
    alpha = 0.05
    min_n = 3

    allelic_freqs1 = segment1['Allelic Frequencies']
    allelic_freqs2 = segment2['Allelic Frequencies']

    # check if there is sufficient number of observations, if no, segments are compared only by copy number
    if len(allelic_freqs1) < min_n and len(allelic_freqs2) < min_n:
        return True

    if len(allelic_freqs1) < min_n or len(allelic_freqs2) < min_n:
        return False

    statistic, p_value = stats.ttest_ind(allelic_freqs1, allelic_freqs2, equal_var=False)

    return p_value > alpha


# count LST score for input sample
def count_lsts(data, LST_SMb=11*Mb, S_small=3*Mb):
    lsts = 0
    for index, segment in data.iterrows():
        # not last segment in profile?
        if index != len(data) - 1:

            _next = data.loc[index+1]
            if segment['Length'] >= LST_SMb and _next['Length'] >= LST_SMb and _next['Chromosome'] == segment['Chromosome'] and _next['Arm'] == segment['Arm'] \
                and _next['Start'] - segment['End'] < S_small:

                lsts += 1

    return lsts
