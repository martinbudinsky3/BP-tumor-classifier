import pandas as pd

def process_centromeres_data(centromeres):
    centromeres = centromeres[centromeres['type'] == 'centromere']
    centromeres = centromeres.loc[:, ['chrom', 'chromStart','chromEnd', 'size']]
    centromeres.rename(columns = { 'chrom': 'Chromosome', 'chromStart': 'Start', 'chromEnd': 'End', 'size': 'Length' }, inplace = True)
    centromeres.iloc[:, 0] = centromeres.loc[:, 'Chromosome'].apply(lambda x: x.replace('chr', ''))
    centromeres = centromeres.reset_index(drop=True)
    
    return centromeres