#!/usr/bin/env python

# In[1]:


import numpy as np
import pandas as pd
import os

data = "/w210/ucb-w210-cap-drug-response/data/"
pfile = "/w210/ucb-w210-cap-drug-response/Pickle_Files/"

all_joined2 = pd.read_csv(data + "model1_to_score.csv")
scoreddf=pd.read_csv(data + "scored_model1.csv")
drug_d=pd.read_csv(pfile + "Screened_Compounds.csv",sep=',')

scored=pd.concat([all_joined2,scoreddf],axis=1)[['COSMIC_DRUG_ID','Prediction']]
scored["DRUG_ID"]=(scored['COSMIC_DRUG_ID'].str.extract(pat = '(["_"].+)'))
scored["DRUG_ID"]=scored["DRUG_ID"].str.replace('_','').astype(int)

scored["COSMIC_ID"]=(scored['COSMIC_DRUG_ID'].str.extract(pat = '(.+["_"])'))
scored["COSMIC_ID"]=scored["COSMIC_ID"].str.replace('_','').astype(int)


tdf=pd.read_csv(pfile + "Drug_Response_threshold.csv")
#tdf=tdf.drop(['Screened Compounds:'],axis=1)
print(tdf.head(2))
scored_resp=pd.merge(scored,tdf,how='left',on='DRUG_ID')
scored_resp=scored_resp.rename(columns={'Prediction':'Predicted_IC50'})
scored_resp['Predicted_Resp']=np.where(scored_resp['Predicted_IC50']>=scored_resp['IC50_Threshold'],"R","S")
scored_resp['Predicted_Resp']=np.where(scored_resp['IC50_Threshold'].isnull(),"UNK",scored_resp['Predicted_Resp'])
print(scored_resp['Screened Compounds:'].head(5))
scored_resp['CONDITION'] = "SCLC"
scored_resp["PATHWAY"] = "RTK signaling"
scored_resp['MODEL'] = 'SCLC'
scored_resp.rename(columns= {'COSMIC_DRUG_ID':'COSMIC_DRUG_ID' \
                            ,'Predicted_IC50':'LN_IC50' \
                            ,'DRUG_ID':'DRUG_ID' \
                            ,'COSMIC_ID':'PATIENT_ID' \
                            ,'IC50_Threshold':'THRESHOLD' \
                            ,'Predicted_Resp':'BINARY_RESPONSE' \
                            ,'Screened Compounds:':'DRUG_NAME' \
                            } \
                    ,inplace=True)
scored_resp.replace('','None', inplace=True)
print(scored_resp['DRUG_NAME'].head(2))
print(scored_resp.columns)
scored_resp.to_csv(data + "complete_score.csv", header = False, na_rep = 'None', index=False)
print(scored_resp[scored_resp['PATIENT_ID'] ==688014].head(2))

print("Finished Post-Processing")
