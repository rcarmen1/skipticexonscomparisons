Negative strand transcription  <--<--<--<--<--<--<--<--

import pandas as pd
import numpy as np
import re
#junction HYOU1 chr11:119,052,112-119,052,711 sj 119,052,636
#junction XPNPEP1 chr10:109,884,115-109,888,097 sj 109,888,048
#junction MYBBP1A chr17:4,545,153-4,545,622 sj 4,545,609
'''gene='HYOU1'
coordinates='chr11:119,052,112-119,052,711'
startj=119052636'''

'''gene='XPNPEP1'
coordinates='chr10:109,884,115-109,888,097'
startj=109888048'''

gene='MYBBP1A'
coordinates='chr17:4,545,153-4,545,622'
startj=4545609

coordinates = re.sub(",", "", coordinates)
print(coordinates)
name=f'5psi_r_{gene}'

url=f'https://snaptron.cs.jhu.edu/gtexv2/snaptron?regions={coordinates}&rfilter=samples_count%3E50'
print(url)
url_df=pd.read_table(url,sep='\t')
#print(url_df.head(2))

subset_url_df = url_df[url_df.end==startj] # .end (reverse junction)
print(subset_url_df.shape)
print(subset_url_df)

# define lambda function to convert string to dictionary
string_to_dict = lambda x: {int(k):int(v) for k,v in [item.split(':') for item in x.split(',')[1:]]} if x.startswith(',') \
else {k:int(v) for k,v in [item.split(':') for item in x.split(',')]} 

# apply lambda function to 'samples' column and assign result to new column 'sample_dict'
subset_url_df['samples'].apply(string_to_dict)
subset_url_df['samples_dict']=subset_url_df['samples'].apply(string_to_dict)

subset_url_df.shape
for n in subset_url_df['samples_dict']:
    print(len(n),type(n))

    
# read the sample TSV file into a DataFrame (sample == df) 
df = pd.read_csv('gtxv2_samples.tsv', sep='\t',low_memory=False)

# print the name of all columns and its types
#print(df.dtypes)
#print(df.shape)
#print(df.head(3))

#n=df['rail_id'].map(subset_url_df['samples'][2])
#printn.head(3)

#map Splice junction with rai_id
i=0
for index, row in subset_url_df.iterrows():
    i=subset_url_df.loc[index,'start']
    colname = f'5junction_{gene}_{i}'
    df[colname]=df['rail_id'].map(subset_url_df.loc[index,'samples_dict'])
    
print(df.head(5))

smts = df['SMTS'].unique().tolist()
print(smts)
subset_df = df[df['SMTS'].isin(['Brain','Nerve','Pituitary'])]
print(subset_df.shape)

df.to_csv(f'{name}_{startj}.csv') 
subset_df.to_csv(f'brain_{name}_{startj}.csv') 
