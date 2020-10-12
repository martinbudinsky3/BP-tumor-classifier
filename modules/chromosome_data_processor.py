import pandas as pd


def process_chromosome_data(data):    
    data = reshape_data(data)
    data = process_values(data)
    data = rename_columns(data)
    data = remove_rows(data)
    data = cast_column_types(data)
    
    return data
    
    
def reshape_data(data):
    data = data.loc[:, 'Chromosome Region':'Length']
    
    data['Start'] = pd.Series('', index=data.index)
    data['End'] = pd.Series('', index=data.index)
    
    return data
    
    
def process_values(data):
    data[['Chromosome Region', 'Start', 'End']] = data['Chromosome Region'].apply(process_chromosome_reg_col)
    
    data['Chromosome Region'] = data['Chromosome Region'].apply(lambda x: x.strip('chr'))
    data['Start'] = data['Start'].apply(lambda x: x.replace(',', ''))
    data['End'] = data['End'].apply(lambda x: x.replace(',', ''))

    category_to_number = {'Event' : {'Big Loss' : 0,
                                'CN Loss' : 1,
                                'CN Gain' : 3,
                                'High Copy Gain' : 4}
                     }
    data.replace(category_to_number, inplace=True)
    
    return data
    
    
def process_chromosome_reg_col(text):
    chromosome, rest = text.split(':')
    
    start, end = rest.split('-')
    return pd.Series([chromosome, start, end])


def rename_columns(data):
    data.rename(columns = {'Chromosome Region':'Chromosome'}, inplace = True)
    data.rename(columns = {'Event':'Copy Number'}, inplace = True)
    
    return data


def remove_rows(data):
    data = data[(data['Copy Number'] != 'LOH') & (data['Copy Number'] != 'Allelic Imbalance') & (data['Copy Number'] != 'Homozygous Copy Loss')]
    
    return data


def cast_column_types(data):
    data['Start'] = data['Start'].astype(str).astype('int64')
    data['End'] = data['End'].astype(str).astype('int64')
    data['Copy Number'] = data['Copy Number'].astype(str).astype('int64')
        
    return data