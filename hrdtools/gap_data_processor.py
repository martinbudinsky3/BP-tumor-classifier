import pandas as pd

class GapDataProcessor:
    """
    Class for preprocessing data of gap regions such as coordinates of telomeres and centromeres 
        
    Methods
    -------
    get_centromeres()
        Returns preprocessed data of centromeres
        
    get_telomeres()
        Returns preprocessed data of telomeres
    """
    
    def __init__(self, filename):
        """
        Parameters
        ----------
        filename : str
            Path to file with gap data
        """
        
        self.data = pd.read_csv(filename, sep='\t')
        self.data = self.process_data(self.data)
        
    
    # method for common preprocessing of all type of gap regions
    def process_data(self, data):
        data = data.loc[(data['type'] == 'centromere') | (data['type'] == 'telomere'), ['chrom', 'chromStart','chromEnd', 'type']]
        data.rename(columns = { 'chrom': 'Chromosome', 'chromStart': 'Start', 'chromEnd': 'End' }, inplace = True)
        data.iloc[:, 0] = data.loc[:, 'Chromosome'].apply(lambda x: x.replace('chr', ''))

        return data
    
    
    def get_centromeres(self):
        """
        Method that returns pandas.DataFrame containing centromeres coordinates
        """
        
        centromeres = self.data.copy()
        centromeres = centromeres.loc[centromeres['type'] == 'centromere']
        centromeres.drop(columns='type', inplace=True)
        centromeres.reset_index(drop=True, inplace=True)

        return centromeres
    
    
    def get_telomeres(self):
        """
        Method that returns pandas.DataFrame containing telomeres coordinates
        """
        
        telomeres = self.data.copy()
        telomeres = telomeres.loc[telomeres['type'] == 'telomere']
        telomeres.drop(columns='type', inplace=True)
        telomeres = telomeres.reset_index(drop=True)

        return telomeres