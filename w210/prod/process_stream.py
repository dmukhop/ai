#!/usr/bin/env python
import json
import os
from pyhive import presto
import numpy as np
import pandas as pd
import os

data = "/w210/ucb-w210-cap-drug-response/data/"
pfile = "/w210/ucb-w210-cap-drug-response/Pickle_Files/"

model_df = pd.read_csv(data + "model1_to_score.csv")
cursor = presto.connect(host = 'presto', port='8080').cursor()
sql_txt = 'SELECT value from raw_input'
print(sql_txt)
s3 = []
cursor.execute(sql_txt)
while True:
    a = cursor.fetchone()
    if a == None :
        break
    start = str(a).find('"COSMIC_ID":' )+13
    end = start + 6
    print(start, end, str(a)[start:end] )
    s3.append(str(a)[start:end] )
cursor.close()
model_df["COSMIC_ID"] =  model_df["COSMIC_ID"]=(model_df['COSMIC_DRUG_ID'].str.extract(pat = '(.+["_"])'))
model_df["COSMIC_ID"]=model_df["COSMIC_ID"].str.replace('_','').astype(int) 
print(model_df["COSMIC_ID"].head(5))
out_model_df = model_df[model_df.COSMIC_ID.isin(s3)]
out_model_df.to_csv(data + "model1_to_score_mini.csv", header = False, na_rep = 'None', index=False)
print('Filtered data for :',s3)


