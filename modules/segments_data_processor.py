import pandas as pd


class SegmentsDataProcessor:
    
    def __init__(self, filename):
        self.data = pd.read_csv(filename, sep='\t', comment='#')
        self.data = self.process_data(self.data)
    
    def process_data(self, data):
        data = self.reshape_data(data)
        data = self.process_values(data)
        data = self.rename_chromosome_column(data)
        data = self.cast_column_types(data)
        data = self.remove_Y_rows(data)
        
        return data
        
        
    def get_cnv_segments(self, sex='female'):    
        cnv_data = self.data.copy()
        cnv_data = self.remove_not_cnv_rows(cnv_data)
        cnv_data = self.rename_event_column(cnv_data)
        cnv_data = self.transform_cn_values(cnv_data)
        cnv_data = self.cast_cn(cnv_data)
        cnv_data.reset_index(inplace=True, drop=True)
        
        return cnv_data
    
    
    def get_ai_segments(self):    
        ai_data = self.data.copy()
        ai_data = self.remove_not_ai_rows(ai_data)
        ai_data = self.drop_event_column(ai_data)
        ai_data.reset_index(inplace=True, drop=True)
        
        return ai_data
    
    
    def get_loh_segments(self):
        loh_data = self.data.copy()
        loh_data = self.remove_not_loh_rows(loh_data)
        loh_data = self.drop_event_column(loh_data)
        loh_data.reset_index(inplace=True, drop=True)
        
        return loh_data


    def reshape_data(self, data):
        data = data.loc[:, 'Chromosome Region':'Length']

        data.loc[:,'Start'] = pd.Series('', index=data.index)
        data.loc[:,'End'] = pd.Series('', index=data.index)

        return data


    def process_values(self, data):
        data[['Chromosome Region', 'Start', 'End']] = data['Chromosome Region'].apply(self.process_chromosome_reg_col)
        
        data.loc[:,'Chromosome Region'] = data.loc[:,'Chromosome Region'].apply(lambda x: x.strip('chr'))
        data.loc[:,'Start'] = data.loc[:,'Start'].apply(lambda x: x.replace(',', ''))
        data.loc[:,'End'] = data.loc[:,'End'].apply(lambda x: x.replace(',', ''))
        data.loc[:,'Length'] = data.loc[:,'Length'].apply(lambda x: x - 1)

        return data


    def process_chromosome_reg_col(self, text):
        chromosome, rest = text.split(':')
        start, end = rest.split('-')
        
        return pd.Series([chromosome, start, end])


    def rename_chromosome_column(self, data):
        data.rename(columns = {'Chromosome Region':'Chromosome'}, inplace = True)

        return data


    def remove_Y_rows(self, data):
        data = data.loc[data['Chromosome'] != 'Y']

        return data


    def cast_column_types(self, data):
        data.loc[:,'Start'] = data.loc[:,'Start'].astype(str).astype('int64')
        data.loc[:,'End'] = data.loc[:,'End'].astype(str).astype('int64')

        return data

    
    def transform_cn_values(self, data):
        category_to_number = {'Copy Number' : {'Big Loss' : 0,
                                'CN Loss' : 1,
                                'CN Gain' : 3,
                                'High Copy Gain' : 4}
                         }
        data.replace(category_to_number, inplace=True)
        
        return data
    
    
    def remove_not_cnv_rows(self, data):
        data = data.loc[(data['Event'] != 'LOH') & (data['Event'] != 'Allelic Imbalance') & (data['Event'] != 'Homozygous Copy Loss')]

        return data
    
    
    def rename_event_column(self, data):
        data.rename(columns = {'Event':'Copy Number'}, inplace = True)

        return data
    
    
    def cast_cn(self, data):
        data.loc[:,'Copy Number'] = data.loc[:,'Copy Number'].astype(str).astype('int64')

        return data
    
    
    def remove_not_ai_rows(self, data):
        data = data.loc[data['Event'] == 'Allelic Imbalance']

        return data
    
    
    def remove_not_loh_rows(self, data):
        data = data.loc[data['Event'] == 'LOH']

        return data
    
    
    def drop_event_column(self, data):
        data.drop(columns='Event', inplace=True)
        
        return data
    
    
class SegmentsDataProcessor2(SegmentsDataProcessor):
    
    def __init__(self, filename, sample_name):
        self.sample_name = sample_name
        self.data =  pd.read_csv(filename, sep='\t', header=None)
        self.data = self.process_data(self.data)
        

    def process_data(self, data):
        data = self.get_sample_data(data)
        data = self.reshape_data_2(data)
        data = self.add_header(data)
        
        return super().process_data(data)


    def get_sample_data(self, data):
        data = data.loc[data[0] == self.sample_name] # 'P6.Rec7'

        return data


    def reshape_data_2(self, data):
        data = data.loc[:, '1':'3']

        data.loc[:,'4'] = pd.Series('', index=data.index)
        data.loc[:,'5'] = pd.Series('', index=data.index)

        return data


    def add_header(self, data):
        data.columns = ['Chromosome Region', 'Event', 'Length', 'Start', 'End']

        return data