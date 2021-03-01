import pandas as pd
from tumor_data_processor import *


def process_tumor_data_with_sex(data, sex='female'):
    data = reshape_data(data)
    data = process_values(data)
    data = rename_columns(data)
    data = remove_rows_with_other_events(data)
    data = cast_column_types(data)
    if sex == 'male':
        data = correct_length(data)
    data.reset_index(inplace=True, drop=True)
    
    return data


def remove_rows_with_other_events(data):
    data = data[(data['Copy Number'] != 'LOH') & (data['Copy Number'] != 'Allelic Imbalance') & (data['Copy Number'] != 'Homozygous Copy Loss')]
    data = data[data['Chromosome'] != 'Y']
    
    return data


def correct_length(data):
    data.loc[data['Chromosome'] == 'X', 'Copy Number'] = data.loc[data['Chromosome'] == 'X', 'Copy Number'].apply(lambda x: x - 1)  
    return data