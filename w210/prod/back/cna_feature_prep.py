#!/usr/bin/env python
# coding: utf-8

# In[1]:
print(1)

import numpy as np
import os,sys
import re
import pandas as pd
from io import StringIO 
#import seaborn as sns
import matplotlib.pyplot as plt
from mapper import expand, parse_mapping_table, apply_mappers
data = "/w210/ucb-w210-cap-drug-response/data/"

# In[2]:
print(2)

GDSC_CNA = data + "Gene_level_CN.csv"

gdsc = pd.read_csv(data + 'Gene_level_CN.csv')
#gdsc.head(3)


# In[3]:
print(3)

gdsc.set_index('gene',inplace = True)
gdsc.drop(["chr","start","stop"],inplace=True,axis=1)
gdsc.columns = gdsc.iloc[0,:]
gdsc = gdsc.iloc[1:,:]
gdsc.columns.name = None
# replace 2001-12-01 with DEC1 and get remove gene names converted to datetimes
gdsc.index.values[37778] = "DEC1"
df_size = gdsc.shape[0]
ndxs=pd.Series(gdsc.index).apply(lambda x : type(x) == str)
gdsc = gdsc.loc[gdsc.index.values[ndxs[ndxs].index],:]
#print(df_size - gdsc.shape[0],"gene IDs excluded due to string to datetime conversion in Excel.")

gdsc.index.name = "gene_id"
ids = gdsc.index
ids = list(set(ids[ids.duplicated()]))
#print("Strings containing duplicated gene IDs:",gdsc.loc[ids,:].shape[0])


# In[4]:
print(4)

def define_avg_ploidy(col):
    n,pl = 0,0
    CN_non_disrupted = []
    for code in col.values:
        if not code == "-1,-1,-,-":
            [max_cn,min_cn,zygosity,disruption] = code.split(",")
            n+=1
            cn = (int(max_cn)+int(min_cn))*0.5
            pl += cn
            if not disruption == "D":
                CN_non_disrupted.append((cn))
    return pd.Series({"avg_pl":pl/n , "median_pl":np.median(CN_non_disrupted)})


# In[5]:
print(5)

GDSC_Ploidies = data + "PICNIC_average_ploidies.tsv"
GDSC_Ploidies = pd.read_csv(GDSC_Ploidies,sep = "\t")
GDSC_Ploidies.drop("#sample_name",axis = 1, inplace= True)
GDSC_Ploidies.set_index("sample_id",inplace=True)
#print(GDSC_Ploidies.shape[0])
GDSC_Ploidies.dropna(inplace=True)
#print(GDSC_Ploidies.shape[0])

est_ploidies = gdsc.apply(define_avg_ploidy).T
df_ploidies = pd.DataFrame.from_dict({"est. avg. ploidy from CN profile":est_ploidies["avg_pl"],"PICNIC avg. pl.":GDSC_Ploidies["average_ploidy"],
                                     "est. median. ploidy":est_ploidies["median_pl"]})


# In[6]:
print(6)

#plt.figure(figsize=(20,5))
#plt.subplot(131)
#tmp = plt.hist(est_ploidies["avg_pl"],bins=30)
#plt.title("est. avg. ploidy from CN profile")
#plt.subplot(132)
#tmp = plt.hist(sorted(list(GDSC_Ploidies["average_ploidy"].values)),bins=30)
#plt.title("PICNIC avg. pl.")
#plt.subplot(133)
#tmp = plt.hist(est_ploidies["median_pl"],bins=30)
#plt.title("est. median ploidy")

#tmp = df_ploidies.plot.scatter(x = "est. avg. ploidy from CN profile",y="PICNIC avg. pl.")


# In[7]:


#est_ploidies.head()


# In[8]:


#df_ploidies.describe()


# In[9]:


# PICNIC average ploidy vs estimated copy-neutral 
#tmp = df_ploidies.boxplot(column="PICNIC avg. pl.", by = "est. median. ploidy" )


# ## Log Transform

# In[10]:
print(7)

estimated_CN = est_ploidies["median_pl"].to_dict()
#estimated_CN['1287381']


# In[11]:


num_marker_thr = 5
# to detect 1 copy gains or losses presenting at CCF >= 0.3
pos_seg_mean_thr = 0.20
neg_seg_mean_thr = -0.23 


def CN2log2R(col, median_ploidy=2 ):
    # this is fr GDSC only
    lRs = []
    genes = col.index.values
    for code in col.values:
        if not code == "-1,-1,-,-":
            [max_cn,min_cn,zygosity,disruption] = code.split(",")
            if int(max_cn) == 0:
                lRs.append(-4.32) # CN=0 with 95% purity
            else:
                max_lR = np.log2(float(max_cn)/median_ploidy)
                if not disruption == "D":
                    lRs.append(max_lR)
                else:
                    if int(min_cn) == 0:
                        min_lR = -4.32
                    else:
                        min_lR = np.log2(float(min_cn)/median_ploidy)
                    if abs(min_lR) > abs(max_lR):
                        lRs.append(min_lR)
                    else:
                        lRs.append(max_lR)
                
        else:
            lRs.append(np.NaN)
    return pd.Series(dict(zip(genes, lRs)))

def clean_logR(logR_value, pos_seg_mean_thr, neg_seg_mean_thr):
    if logR_value >= pos_seg_mean_thr:
        return logR_value 
    elif logR_value <= neg_seg_mean_thr:
            return logR_value 
    else:
        return 0
    
def handle_dups(df,corr_thr = 0.75):
    '''Detect dupliated row IDs. Merge 2 or more rows with the same ID, 
    if averaged correlation in all pairvise comparision is >= corr_thhr;\n
    otherwise drop all duplicates.  Keeps abs. max value (negative preferred).'''
    dups = df.index
    dups = list(set(dups[dups.duplicated()]))
    if len(dups)==0:
        print("No duplicated row IDs. Do nothing.")
        return df
#    print(len(dups), "duplicated IDs in",df.loc[dups,:].shape[0],"rows found.")
    dups_merge = [] # if corr > corr_thr
    dups_remove = [] # corr < 
    for dup in dups:
        r = df.loc[dup,:].T.corr()
        n_dups = df.loc[dup,:].shape[0]
        r_avg = []
        for i in range(0,n_dups):
            for j in range(i+1,n_dups):
                r_avg.append(r.iloc[i,j])
        if np.average(r_avg) < corr_thr :
            #print(dup,r_avg, n_dups)
            dups_remove.append(dup)
        else:
            dups_merge.append(dup)
    
    # remove not similar duplicates
    df_size = df.shape[0]
    df = df.loc[~df.index.isin(dups_remove),:]
#    print("duplicate rows removed due to low correlation of duplicated profiles",df_size -df.shape[0] )
    df_size = df.shape[0]
    
    # merge simialr duplicates
    d1 = df.loc[~df.index.isin(dups_merge),:]
    d2 = df.loc[dups_merge,:]
    d2 = d2.groupby(d2.index).agg(lambda x: -max(-x.max(),-x.min(),key= abs))
    df = pd.concat([d1,d2])
    df.sort_index(inplace=True)
#    print("Merged ",df_size-df.shape[0]+len(dups_merge),"duplicated rows into",len(dups_merge),"rows")
    return df


# In[12]:
print(8)

gdsc = gdsc.apply(lambda x : CN2log2R(x,estimated_CN[x.name] ))
# drop genes without any determined value
gdsc = gdsc.dropna(axis=0,how="all")
# fill with zeroes the remaining ones
gdsc.fillna(0,inplace=True)
#gdsc.head(3)


# In[13]:
print(9)

gdsc = gdsc.applymap(lambda x :  clean_logR(x, pos_seg_mean_thr, neg_seg_mean_thr))


# In[14]:


NCBI = pd.read_csv(data + "Homo_sapiens.gene_info",sep = "\t")
NCBI = NCBI[["#tax_id","GeneID","Symbol","Synonyms","type_of_gene"]]
NCBI = NCBI.loc[NCBI["#tax_id"] == 9606]
NCBI = NCBI.loc[NCBI["type_of_gene"] != "unknown"]
ncbi_symbols = parse_mapping_table(NCBI, "Symbol","GeneID")


# In[15]:
print(10)

ncbi_synonyms = expand(NCBI[["Synonyms","GeneID"]],column="Synonyms",sep="|") 
ncbi_synonyms = parse_mapping_table(ncbi_synonyms, "Synonyms","GeneID")


# In[16]:


#NCBI.head()


# In[17]:


gdsc,query2target,not_mapped  = apply_mappers(gdsc, ncbi_symbols, ncbi_synonyms, verbose = True,handle_duplicates = "keep")
#gdsc.head(3)


# In[18]:


gdsc = handle_dups(gdsc,corr_thr = 0.75)


# In[19]:


gdsc.index.name = "gene_id"
gdsc = gdsc.T.sort_index().T
#gdsc.head()


# In[20]:
print(11)

gdsc.to_csv(data + "GDSC.Segment_Mean.CNA.tsv", sep = "\t",header=True,index=True)


# In[21]:


cna_segment_mean = pd.read_csv(data + 'GDSC.Segment_Mean.CNA.tsv', index_col='gene_id', sep= "\t")


# In[22]:


#cna_segment_mean.head()


# In[23]:
print(12)

cna_segment_mean = cna_segment_mean.transpose()
cna_segment_mean.index.names = ['COSMIC_ID']
#cna_segment_mean.head()


# In[24]:


cna_segment_mean.to_csv(data + "cna_scored_transposed.tsv",sep = "\t")

print("End of CNA Processing")
# In[ ]:




