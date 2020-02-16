#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import os

data = "/w210/ucb-w210-cap-drug-response/data/"
pfile = "/w210/ucb-w210-cap-drug-response/Pickle_Files/"

all_joined2 = pd.read_csv(data + "model1_to_score.csv")
scoreddf=pd.read_csv(data + "scored_model1.csv")

scored=pd.concat([all_joined2,scoreddf],axis=1)[['COSMIC_DRUG_ID','Prediction']]
scored["DRUG_ID"]=(scored['COSMIC_DRUG_ID'].str.extract(pat = '(["_"].+)'))
scored["DRUG_ID"]=scored["DRUG_ID"].str.replace('_','').astype(int)

scored["COSMIC_ID"]=(scored['COSMIC_DRUG_ID'].str.extract(pat = '(.+["_"])'))
scored["COSMIC_ID"]=scored["COSMIC_ID"].str.replace('_','').astype(int)


tdf=pd.read_csv(pfile + "Drug_Response_threshold.csv")
tdf=tdf.drop(['Screened Compounds:'],axis=1)

scored_resp=pd.merge(scored,tdf,how='left',on='DRUG_ID')
scored_resp=scored_resp.rename(columns={'Prediction':'Predicted_IC50'})
scored_resp['Predicted_Resp']=np.where(scored_resp['Predicted_IC50']>=scored_resp['IC50_Threshold'],"S","R")
scored_resp['Predicted_Resp']=np.where(scored_resp['IC50_Threshold'].isnull(),"UNK",scored_resp['Predicted_Resp'])

scored_resp.to_csv(data + "complete_score.csv")

print("Finished Post-Processing")
