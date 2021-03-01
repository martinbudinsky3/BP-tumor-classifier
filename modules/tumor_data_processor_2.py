import pandas as pd
from tumor_data_processor import *


def process_tumor_data_2(data, sex='female'):
    data = remove_other_tumors_data(data)
    data = reshape_data_2(data)
    data = add_header(data)
    data = process_values(data)
    data = rename_columns(data)
    data = remove_rows_with_other_events_2(data)
    data = cast_column_types(data)
    if sex == 'male':
        data = correct_length_2(data)
    data.reset_index(inplace=True, drop=True)
    
    return data


def remove_other_tumors_data(data):
    data = data[data[0] == 'P6.Rec7']
    
    return data


def reshape_data_2(data):
    data = data.loc[:, '1':'3']
    
    data['4'] = pd.Series('', index=data.index)
    data['5'] = pd.Series('', index=data.index)
    
    return data


def add_header(data):
    data.columns = ['Chromosome Region', 'Event', 'Length', 'Start', 'End']
    
    return data


def remove_rows_with_other_events_2(data):
    data = data[(data['Copy Number'] != 'LOH') & (data['Copy Number'] != 'Allelic Imbalance') & (data['Copy Number'] != 'Homozygous Copy Loss')]
    data = data[data['Chromosome'] != 'Y']
    
    return data


def correct_length_2(data):
    data.loc[data['Chromosome'] == 'X', 'Copy Number'] = data.loc[data['Chromosome'] == 'X', 'Copy Number'].apply(lambda x: x - 1)  
    
    return data