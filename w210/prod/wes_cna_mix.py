#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
#import seaborn as sns
import matplotlib.pyplot as plt
from itertools import islice
import os
from sklearn.decomposition import PCA
from sklearn.externals import joblib
# !cat /proc/sys/vm/overcommit_memory
# !echo 1 > /proc/sys/vm/overcommit_memory

data = "/w210/ucb-w210-cap-drug-response/data/"
pfile = "/w210/ucb-w210-cap-drug-response/Pickle_Files/"

# In[2]:


#Read Genomic data


# In[3]:


gen_data=pd.read_csv(pfile + "genomic_features.tsv",sep="\t",nrows=10)

#Optimize datatypes to save memory
gene_dict=gen_data.dtypes.apply(lambda x: x.name).to_dict()
new_gene_dict=gene_dict.copy()
vlist=['COSMIC_ID','Sample Name','TISSUE_FACTOR']
i=1
for k in new_gene_dict.keys():
        if k in vlist:
            i=0
        else:
            new_gene_dict[k]='int16'
#         print(k,new_wes_dict[k])

genomic_features = pd.read_csv(pfile + "genomic_features.tsv",sep = "\t",dtype=new_gene_dict)

#gen_data.head(2)


# In[4]:


# Cell line details


# In[5]:


celld=pd.read_csv(data + 'Cell_Lines_Details.csv')
#celld.head(2)


# In[6]:


# Read in drug data
drug_d=pd.read_csv(pfile + "Screened_Compounds.csv",sep=',')
#drug_d.head()


# In[7]:


# Read in drug response data
dose_response_df_all = pd.read_csv(data + 'v17.3_fitted_dose_response_binary.csv', index_col='IC50_RESULTS_ID')
dose_response_df=pd.merge(dose_response_df_all,drug_d,how='left',on='DRUG_ID')
dose_response_df['TARGET_PATHWAY']=dose_response_df['TARGET_PATHWAY'].fillna("UNK")
#dose_response_df.head(2)


# In[8]:


#Read in cancer labels, clean data
tcga_labels=pd.read_csv(pfile + "TCGA_Labels.csv")
celld=pd.read_csv(data + 'Cell_Lines_Details.csv')
celld=celld[['COSMIC identifier','GDSC\nTissue descriptor 1','GDSC\nTissue\ndescriptor 2','Cancer Type\n(matching TCGA label)']]
celld_desc=pd.merge(celld,tcga_labels,how='outer',left_on='Cancer Type\n(matching TCGA label)',right_on="Sigle",indicator=True)
celld_desc['Program']=celld_desc['Program'].fillna('UNK')
celld_desc=celld_desc.rename(columns={'COSMIC identifier':'COSMIC_ID','GDSC\nTissue descriptor 1':'GDSC_DESC1','GDSC\nTissue\ndescriptor 2':'GDSC_DESC2','Cancer Type\n(matching TCGA label)':'TCGA_Label', 'Program':'Cancer_Type'                                     })
dose_resp_caname=pd.merge(dose_response_df,celld_desc,how='left',on='COSMIC_ID')
dose_response_trimmed_df = dose_resp_caname[['COSMIC_ID', 'DRUG_ID', 'LN_IC50', 'BINARY_RESPONSE','TCGA_Label','Cancer_Type','TARGET_PATHWAY']]
#dose_response_trimmed_df.head(2)


# In[9]:


# WES prep


# In[10]:


# Read Wes data first few rows
wes_scored_temp = pd.read_csv(data + "wes_binary_scored_transposed.tsv",sep = "\t",nrows=10)
#wes_scored_temp.head()


# In[11]:


# Read full WES data with right data types
wes_dict=wes_scored_temp.dtypes.apply(lambda x: x.name).to_dict()
new_wes_dict=wes_dict.copy()
i=1
for k in new_wes_dict.keys():
        if k=='COSMIC_ID':
            i=0
        else:
            new_wes_dict[k]='int16'
#         print(k,new_wes_dict[k])

wes_scored = pd.read_csv(data + "wes_binary_scored_transposed.tsv",sep = "\t",dtype=new_wes_dict)


# In[12]:


wes_scored=wes_scored.rename(columns={"Unnamed: 0":"COSMIC_ID"})


# In[13]:


#Subset one pathway RTK signaling and lung cancer - THIS MAY NEED TO BE PARAMETERS
rtk=dose_resp_caname[(dose_resp_caname['TARGET_PATHWAY']=='RTK signaling')  & (dose_resp_caname['GDSC_DESC1']=='lung_SCLC')]
# rtk['Cancer_Type'].value_counts()
#rtk.shape


# In[14]:


rtk_wes_joined_df = pd.merge(wes_scored, rtk,how='inner', on='COSMIC_ID')
print(rtk_wes_joined_df.shape)
print(rtk.shape)


# In[15]:


rtk_wes_joined_df['COSMIC_DRUG_ID']=rtk_wes_joined_df['COSMIC_ID'].map(str)+"_"+rtk_wes_joined_df['DRUG_ID'].map(str)
wes_subset=rtk_wes_joined_df[['COSMIC_DRUG_ID','347733','440560','4633','2019','5706','192683','1063','2220','4829','5137','7767','7982','9994','23033','51585','55051','57509','65083','84695','119395','148137','387266']]


# In[16]:


#CNA


# In[17]:


cna_scored_temp = pd.read_csv(data + "cna_scored_transposed.tsv",sep = "\t", nrows=10)
#cna_scored_temp.head()


# In[18]:


cna_dict=cna_scored_temp.dtypes.apply(lambda x: x.name).to_dict()
new_cna_dict=cna_dict.copy()
i=1
for k in new_cna_dict.keys():
        if k=='COSMIC_ID':
            i=0
        else:
            new_cna_dict[k]='float16'
#         print(k,new_cna_dict[k])


# In[19]:


cna_scored = pd.read_csv(data + "cna_scored_transposed.tsv",sep = "\t", dtype=new_cna_dict)
#cna_scored.head()
print("PRINT CNA SHAPE:", cna_scored.shape)


# In[20]:


rtk=rtk.drop(['_merge'],axis=1)
rtk_cna_joined_df = pd.merge(cna_scored, rtk, how='inner', on='COSMIC_ID',indicator=True)
rtk_cna_joined_df['COSMIC_DRUG_ID']=rtk_cna_joined_df['COSMIC_ID'].map(str)+"_"+rtk_cna_joined_df['DRUG_ID'].map(str)

#Make a copy with original labels and other fields
rtk_cna_joined_df_copy_with_label=rtk_cna_joined_df.copy()

print("PRINT RTK:", rtk_cna_joined_df.shape)

rtk_cna_joined_df=rtk_cna_joined_df.drop({"LN_IC50","COSMIC_ID","CELL_LINE_NAME","DRUG_ID","_merge","Sigle","GDSC_DESC1","GDSC_DESC2","TCGA_Label","Cancer_Type",'MAX_CONC_MICROMOLAR','MIN_CONC_MICROMOLAR','AUC','RMSE','Z_SCORE','DATASET_VERSION','PUTATIVE_TARGET','DRUG_NAME_x','DRUG_NAME_y','SYNONYMS','TARGET','TARGET_PATHWAY','BINARY_RESPONSE'},axis=1)
rtk_cna_joined_df=rtk_cna_joined_df.set_index("COSMIC_DRUG_ID")


# In[21]:


print("PRINT rtk_cna_joined_df:",rtk_cna_joined_df.shape)


# In[22]:


#Apply PCA on CNA test

# Load from file
joblib_file = pfile + "model1_cna_pca.pkl"
joblib_pca = joblib.load(joblib_file)

# Transform test file
devcnapca=joblib_pca.transform(rtk_cna_joined_df)  

i=1
varlist=[]
for i in range(1, 41):
    var="cna_princ"+str(i)
    varlist.append(var)
    
cna_subset=pd.DataFrame(data=devcnapca,columns=varlist,index=rtk_cna_joined_df.index)
cna_subset=cna_subset.reset_index()
#cna_subset.head()


# In[23]:


#Join the WES and CNA, split COSMIC_DRUG_ID to COSMIC_ID and DRUG_ID
cna_wes_joined=pd.merge(cna_subset,wes_subset,how='inner',on='COSMIC_DRUG_ID')

cna_wes_joined["DRUG_ID"]=(cna_wes_joined['COSMIC_DRUG_ID'].str.extract(pat = '(["_"].+)'))
cna_wes_joined["DRUG_ID"]=cna_wes_joined["DRUG_ID"].str.replace('_','').astype(int)

cna_wes_joined["COSMIC_ID"]=(cna_wes_joined['COSMIC_DRUG_ID'].str.extract(pat = '(.+["_"])'))
cna_wes_joined["COSMIC_ID"]=cna_wes_joined["COSMIC_ID"].str.replace('_','').astype(int)

#cna_wes_joined.head()


# In[24]:


#Transform drug_id to indicator columns
cna_wes_joined = pd.concat([cna_wes_joined, pd.get_dummies(cna_wes_joined['DRUG_ID'], prefix="drug_id_" )],axis=1)


# In[25]:


#Import number of targets with drug id
drug_target=pd.read_csv(pfile + "drugs_num_target.csv")
#drug_target.head(3)


# In[26]:


#Import drugs with target indicator
drug_with_target=pd.read_csv(pfile + "drugs_with_target.csv")
drug_with_target=drug_with_target.drop(['num_targets'],axis=1)
#drug_with_target.head(3)


# In[27]:


#Join all 
cna_wes_gen_joined=pd.merge(cna_wes_joined,genomic_features,how='inner',on='COSMIC_ID')
all_joined=pd.merge(cna_wes_gen_joined,drug_target,how='inner',on='DRUG_ID')
all_joined1=pd.merge(all_joined,drug_with_target,how='inner',on='DRUG_ID')
all_joined2=all_joined1.drop(['Sample Name','TISSUE_FACTOR','DRUG_ID','COSMIC_ID'],axis=1)
#all_joined2.head()


# In[28]:

print("PRINT all_joined2:",all_joined2.shape)
all_joined2.to_csv(data + "model1_to_score.csv")

print("Model Preprocess Complete")

