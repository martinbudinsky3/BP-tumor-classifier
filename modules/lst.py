import pandas as pd
import numpy as np
import scipy.stats as stats

from utils import init_lengths, init_centromeres, init_chromosome_names, Mb

LST_SMbs = [x for x in range(3, 12)]
lengths = init_lengths()
centromeres = init_centromeres()
chromosome_names = init_chromosome_names()


def lst(data, vcf_reader=None, sample_name=None, LST_SMb_param=11):
    data = fill_segments(data)
    dna_index = count_dna_index(data)
        
    # if sample has vcf data coerce by copy numbers and allelic frequencies else only by copy numbers
    if not vcf_reader is None:
        data = count_allele_freqs(data, vcf_reader, sample_name)
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


def insert_row(df, _chr, cn, length, start, end, index):
    normal_segment = pd.DataFrame({
        'Chromosome': [ _chr ],
        'Copy Number': [ cn ],
        'Length': [length],
        'Start': [ start ],
        'End': [ end ]
    })

    return pd.concat([df.iloc[:index], normal_segment, df.iloc[index:]]).reset_index(drop=True)


def fill_segments(data):
    df = data.copy()
    index_df = 0
    normal_cn = 2
    
    for index, row in data.iterrows():

        # first cnv region in chromosome
        if index == 0 or data.loc[index-1, 'Chromosome'] != data.loc[index, 'Chromosome']:
            if row['Start'] != 0:
                df = insert_row(df, row['Chromosome'], normal_cn, row['Start'], 0, row['Start'], index_df)
                index_df += 1

        # not first cnv region in chromosome 
        elif data.loc[index-1, 'End'] != data.loc[index, 'Start']:
            prev = data.loc[index-1]

            df = insert_row(df, row['Chromosome'], normal_cn, row['Start'] - prev['End'], prev['End'], row['Start'], index_df)
            index_df += 1

        # last cnv region in chromosome
        if index == len(data) - 1 or data.loc[index+1, 'Chromosome'] != data.loc[index, 'Chromosome']:

            chr_len = lengths.loc[row['Chromosome'], 'Length']
            if row['End'] != chr_len:      
                df = insert_row(df, row['Chromosome'], normal_cn, chr_len - row['End'], row['End'], chr_len, index_df+1)
                index_df += 1

        index_df += 1

    df = remove_centromeres(df)

    return df


def remove_centromeres(data):
    for _chr in chromosome_names:
        centromere_start = centromeres.loc[centromeres['Chromosome'] == _chr, 'Start'].iloc[0]
        centromere_end = centromeres.loc[centromeres['Chromosome'] == _chr, 'End'].iloc[0]

        # drop segments that start and end in centromere  
        data = data.drop(data[(data['Chromosome'] == _chr) & (data['Start'] >= centromere_start) & (data['End'] <= centromere_end)].index)

        # cut segment that overlaps centromere from both sides
        row_df = data.loc[(data['Chromosome'] == _chr) & (data['Start'] < centromere_start) & (data['End'] > centromere_end)]
        if not row_df.empty:
            row = row_df.iloc[0]
            data.loc[row.name, 'End'] = centromere_start
            data.loc[row.name, 'Length'] = centromere_start - row['Start']        
            data = insert_row(data, _chr, row['Copy Number'], row['End'] - centromere_end, centromere_end, row['End'], row.name + 1)

        # cut segment that ends in centromere  
        data.loc[(data['Chromosome'] == _chr) & (data['End'] >= centromere_start) & (data['End'] <= centromere_end), 'End'] = centromere_start

        # cut segment that starts in centromere  
        data.loc[(data['Chromosome'] == _chr) & (data['Start'] >= centromere_start) & (data['Start'] <= centromere_end), 'Start'] = centromere_end

    data = name_chr_arms(data)

    return data


def name_chr_arms(data):
    data['Arm'] = 'p'
    for _chr in chromosome_names:
        centromere_end = centromeres.loc[centromeres['Chromosome'] == _chr, 'End'].iloc[0]

        data.loc[(data['Chromosome'] == _chr) & (data['Start'] >= centromere_end), 'Arm'] = 'q'

    return data


def count_dna_index(data):
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


# function that checks if vcf record has required quality
def has_quality(record, qual_threshold = 50):
    info = record.INFO
    return record.QUAL > qual_threshold and ('QD' not in info.keys() or info['QD'] > 10.0) and ('MQ' not in info.keys() or info['MQ'] > 40.0) \
        and ('FS' not in info.keys() or info['FS'] < 30.0 ) and ('SOR' not in info.keys() or info['SOR'] < 3.0) \
        and ('MQRankSum' not in info.keys() or info['MQRankSum'] > -12.5) and ('ReadPosRankSum' not in info.keys() or info['ReadPosRankSum'] > -8.0)


# function for counting allelic frequencies for each segment
def count_allele_freqs(data, vcf_reader, sample, qual_threshold = 50):
    data_with_af = data.copy()
    data_with_af['Allele Frequencies'] = [list() for x in range(len(data_with_af.index))]

    for index, row in data_with_af.iterrows():

        try:
            segment_records = vcf_reader.fetch(row['Chromosome'], row['Start'], row['End'])
            allele_freqs = []
        
            for record in segment_records:
                sample_data = record.genotype(sample).data

                if has_quality(record) and sample_data.GT != './.' and sample_data.GT != '0/0' and sample_data.AD != './.' and sample_data.AD != None \
                    and sample_data.AD[1] != 0:

                    allele_freq = sample_data.AD[1] / (sample_data.AD[0] + sample_data.AD[1])
                    allele_freqs.append(allele_freq)

            data_with_af.at[index, 'Allele Frequencies'] = allele_freqs

        except (ValueError, AttributeError) as e:
            continue

    return data_with_af


# insert new segment that was created by linking - segment without allelic frequencies
def insert_row_without_af(df, _chr, cn, length, start, end, arm, index):
    normal_segment = pd.DataFrame({
        'Chromosome': [ _chr ],
        'Copy Number': [ cn ],
        'Length': [length],
        'Start': [ start ],
        'End': [ end ],
        'Arm': [arm]
    })

    return pd.concat([df.iloc[:index], normal_segment, df.iloc[index:]]).reset_index(drop=True)


# insert new segment that was created by linking 
def insert_row_with_af(df, _chr, cn, length, start, end, allele_freqs, arm, index):
    normal_segment = pd.DataFrame({
        'Chromosome': [ _chr ],
        'Copy Number': [ cn ],
        'Length': [length],
        'Start': [ start ],
        'End': [ end ],
        'Allele Frequencies': [allele_freqs],
        'Arm': [arm]
    })

    return pd.concat([df.iloc[:index], normal_segment, df.iloc[index:]]).reset_index(drop=True)


# linking adjacent segments of small filtered out segment - segments without allelic frequencies 
def link_segments_without_af(df, prev, _next):
    df = df.drop(index=prev.name)
    df = df.drop(index=_next.name)

    df = insert_row_without_af(df, prev['Chromosome'], prev['Copy Number'],  _next['End'] - prev['Start'],  prev['Start'], _next['End'], \
                            prev['Arm'], prev.name)

    return df


# linking adjacent segments of small filtered out segment 
def link_segments(df, prev, _next, small):
    df = df.drop(index=prev.name)
    df = df.drop(index=_next.name)

    new_allele_freqs = []
    new_allele_freqs.extend(prev['Allele Frequencies'])
    new_allele_freqs.extend(_next['Allele Frequencies'])
    new_allele_freqs.extend(small['Allele Frequencies'])

    df = insert_row_with_af(df, prev['Chromosome'], prev['Copy Number'],  _next['End'] - prev['Start'],  prev['Start'], _next['End'], new_allele_freqs, \
                            prev['Arm'], prev.name)

    return df


# statistical tests for equality of allelic frequencies of two segments
def have_equal_allele_freqs(segment1, segment2):
    alpha = 0.05
    min_n = 3

    allele_freqs1 = segment1['Allele Frequencies']
    allele_freqs2 = segment2['Allele Frequencies']

    # check if there is sufficient number of observations, if no, segments are compared only by copy number
    if len(allele_freqs1) < min_n and len(allele_freqs2) < min_n:
        return True

    if len(allele_freqs1) < min_n or len(allele_freqs2) < min_n:
        return False

    statistic, p_value = stats.ttest_ind(allele_freqs1, allele_freqs2, equal_var=False)

    if p_value > alpha:
        return True

    return False


# coercing function
def coercing(data, count_af=True, S_small=3*Mb):
    df = data.copy()

    while len(df) > 0:

        # get smallest segment
        row = df[df['Length'] == df['Length'].min()].iloc[0]
        index = row.name

        # filter out?
        if row['Length'] < S_small:

            # not first or last segment of profile?
            if index != 0 and index != len(df) - 1:
                prev = df.loc[ index-1 ]
                _next = df.loc[ index+1 ]

                # can link?
                if prev['Chromosome'] == _next['Chromosome'] and prev['Arm'] == _next['Arm'] and prev['Copy Number'] == _next['Copy Number']:
                    
                    # if sample has vcf data check allelic frequencies else join only based on copy number
                    if not count_af:
                        df = link_segments_without_af(df, prev, _next)
                    elif have_equal_allele_freqs(prev, _next):
                        df = link_segments(df, prev, _next, row)

            # delete small segment
            df = df.drop(index=index).reset_index(drop=True)

        # if there are no small segments left -> end
        else:
            break
        
    return df


def count_lsts(data, LST_SMb=11*Mb, S_small=3*Mb):

    lsts = 0
    for index, row in data.iterrows():

        # not last segment in profile?
        if index != len(data) - 1:

            _next = data.loc[index+1]
            if row['Length'] >= LST_SMb and _next['Length'] >= LST_SMb and _next['Chromosome'] == row['Chromosome'] and _next['Arm'] == row['Arm'] \
                and _next['Start'] - row['End'] < S_small:

                lsts += 1

    return lsts
