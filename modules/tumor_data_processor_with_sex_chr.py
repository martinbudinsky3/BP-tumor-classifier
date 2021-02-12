import pandas as pd
from tumor_data_processor import *


def process_tumor_data_with_sex(data):
    data = reshape_data(data)
    data = process_values(data)
    data = rename_columns(data)
    data = remove_rows_with_other_events(data)
    data = cast_column_types(data)
    data = correct_length(data)
    data.reset_index(inplace=True, drop=True)
    
    return data


def remove_rows_with_other_events(data):
    data = data[(data['Copy Number'] != 'LOH') & (data['Copy Number'] != 'Allelic Imbalance') & (data['Copy Number'] != 'Homozygous Copy Loss')]
    
    return data