#!/usr/bin/env python
# coding: utf-8

# In[1]:
print(1)

import numpy as np
import sys
import re
import pandas as pd
#import seaborn as sns
import matplotlib.pyplot as plt
from mapper import expand, parse_mapping_table, apply_mappers

verbose=False

data = "/w210/ucb-w210-cap-drug-response/data/"
# # Section 1  (Refenece data setup)
# ## Load Onco Gene Data (Reference Information)
# ### The data below is a curated set of genes known to transform a cell into a tumor cell, or a gene known to suppress

# #### OncoKB (Very popular source for cancer gene annotation)

# In[2]:
print(2)
# from http://oncokb.org/#/cancerGenes
OncoKB_mut = pd.read_csv(data+"allAnnotatedVariants.txt",sep ="\t", encoding='latin-1')
OncoKB_gtype = pd.read_csv(data+"allCuratedGenes.txt",sep ="\t")

#OncoKB_gtype.head()

# #### Split out genes into category buckets

# In[3]:
print(3)

og_genes = []
tsg_genes = []
og_tsg_genes = []

for index, row in OncoKB_gtype.iterrows():
    if row['Is Oncogene'] == 'Yes':
        og_genes.append(row['Hugo Symbol'])

    if row['Is Tumor Suppressor Gene'] == 'Yes':
        tsg_genes.append(row['Hugo Symbol'])

    if row['Is Oncogene'] == 'Yes' and row['Is Tumor Suppressor Gene'] == 'Yes':
        og_tsg_genes.append(row['Hugo Symbol'])

OncoKB = {"OG": set(og_genes), "TSG": set(tsg_genes), "OG_TSG": set(og_tsg_genes)}

if verbose:
    print("OncoKB OG:", len(OncoKB["OG"]), "TSG:", len(OncoKB["TSG"]), "OG_TSG:", len(OncoKB["OG_TSG"]))

# ### Bailye (Another source of cancer genes)

# bailye = {"OG":[],"TSG":[],"OG_TSG":[]}
# bailye_drivers = pd.read_excel("../data/mmc1.xlsx", sheet_name="Table S1",
#                                skiprows = 3)
# bailye_drivers = bailye_drivers[["Gene","Tumor suppressor or oncogene prediction (by 20/20+)"]]
# #print(bailye_drivers.shape[0])
# bailye_drivers.drop_duplicates(inplace=True)
# bailye_drivers.columns = ["gene", "role"]
# #print(bailye_drivers.shape[0])
# bailye_drivers.fillna("NA", inplace = True)
# bailye_drivers["role"] = bailye_drivers["role"].apply(lambda x: x.replace("possible ",""))
# bailye_drivers.drop_duplicates(inplace=True)
# #print(bailye_drivers.shape[0])
# bailye_drivers.loc[bailye_drivers.duplicated("gene", keep=False),:]
# #print(len(set(bailye_drivers["gene"].values)))
# bailye["OG"] = bailye_drivers.loc[bailye_drivers["role"] == "oncogene",:]
# bailye["OG"] = set(bailye["OG"]["gene"].values)
# bailye["TSG"] = bailye_drivers.loc[bailye_drivers["role"] == "tsg",:]
# bailye["TSG"] = set(bailye["TSG"]["gene"].values)
# bailye["OG_TSG"] = bailye["OG"].intersection(bailye["TSG"])
# print("Bailey_2018 OG:",len(bailye["OG"]),"TSG:",len(bailye["TSG"]),"both:", len(bailye["OG_TSG"]))

# In[15]:
print(4)

bailye = {"OG": [], "TSG": [], "OG_TSG": []}
bailye_drivers = pd.read_csv(data+ "mmc1_tab.csv")
bailye_drivers = bailye_drivers[["Gene", "Tumor suppressor or oncogene prediction (by 20/20+)"]]
# print(bailye_drivers.shape[0])
bailye_drivers.drop_duplicates(inplace=True)
bailye_drivers.columns = ["gene", "role"]
# print(bailye_drivers.shape[0])
bailye_drivers.fillna("NA", inplace=True)
bailye_drivers["role"] = bailye_drivers["role"].apply(lambda x: x.replace("possible ", ""))
bailye_drivers.drop_duplicates(inplace=True)
# print(bailye_drivers.shape[0])
bailye_drivers.loc[bailye_drivers.duplicated("gene", keep=False), :]
# print(len(set(bailye_drivers["gene"].values)))
bailye["OG"] = bailye_drivers.loc[bailye_drivers["role"] == "oncogene", :]
bailye["OG"] = set(bailye["OG"]["gene"].values)
bailye["TSG"] = bailye_drivers.loc[bailye_drivers["role"] == "tsg", :]
bailye["TSG"] = set(bailye["TSG"]["gene"].values)
bailye["OG_TSG"] = bailye["OG"].intersection(bailye["TSG"])
if verbose:
    print("Bailey_2018 OG:", len(bailye["OG"]), "TSG:", len(bailye["TSG"]), "both:", len(bailye["OG_TSG"]))

# ## NCBI
# Below loads another set of standard reference genes for the human genome. These are all the genes, independent on if they are releated to cancer

# In[5]:
print(5)

NCBI = pd.read_csv(data + "Homo_sapiens.gene_info", sep="\t")
NCBI = NCBI[["#tax_id", "GeneID", "Symbol", "Synonyms", "type_of_gene"]]
NCBI = NCBI.loc[NCBI["#tax_id"] == 9606]
NCBI = NCBI.loc[NCBI["type_of_gene"] != "unknown"]
ncbi_symbols = parse_mapping_table(NCBI, "Symbol", "GeneID")


# ## Helper Methods
#
# Methods below are used to score the impact that the mutations / variants from the WES file have on an individual gene. This allows us to represent all the mutations into a value that can be compared across cell lines

# In[11]:
print(6)

# gene_id - Entrez Gene ID
# df - all mutations in sample
#    columns: "aa_pos","function", "Entrez_gene_ID"
# hotspots - table of known hotsport  "Entrez_gene_ID" and "aa_pos"

def distance_to_hotspot(gene_id, df_gene, hotspots):
    '''Returns minimal distance between point mutation in gene and known hotspot'''
    # consider only GoF because LoF affect the whole protein and not a certian point
    df_gene = df_gene.loc[df_gene["function"] == "point_mutation"]
    #print("DISTANCE_TO_HOTSPOT:",gene_id)
    #print("HOTSPOT:", hotspots.columns )
    df_hotspots = hotspots.loc[hotspots["Entrez_gene_ID"] == gene_id, :]
    dist = 100000
    for s, e in set(df_gene["aa_pos"].values):
        for hs, he in set(df_hotspots["aa_pos"].values):
            dist = min(dist, abs(s - hs), abs(s - he), abs(e - hs), abs(e - he))
    return dist


def mutations2score(df, TSG, OG, hotspots, w=0.05, verbose=False):
    # define gene type
    gene_scores = {}
    for gene_id in list(set(df["Entrez_gene_ID"].values)):
        df_gene = df.loc[df["Entrez_gene_ID"] == gene_id, :]
        # functions of mutations
        mut_functions = set(df_gene["function"].values)
        if gene_id in TSG:
            # check for truncating mutations:
            # score is -1 if any
            if "truncating_mutation" in mut_functions:
                gene_scores[gene_id] = -1
            else:
                # in mutation is not truncating,
                # get the distance to closets hotspot if any
                #print("CHECKING PARAMETERS : Gene_id",gene_id, "df_gene:", df_gene, "hotspot:",hotspots)
                dist = distance_to_hotspot(gene_id, df_gene, hotspots)
                # set abs. maximal of 1/(dist+1) and 2w which corresponds LoF in unannotatedd gene
                gene_scores[gene_id] = -max(1.0 / (dist + 1), 2 * w)
        elif gene_id in OG:
            if "point_mutation" in mut_functions:
                #print("CHECKING PARAMETERS : Gene_id",gene_id, "df_gene:", df_gene, "hotspot:",hotspots)
                dist = distance_to_hotspot(gene_id, df_gene, hotspots)
                if "truncating_mutation" in mut_functions:
                    if verbose:
                        print("GoF at dist", dist, "from hotspot in Oncogene with GoF", file=sys.stderr)
                        print(df_gene, file=sys.stderr)
                    if 1.0 / (dist + 1) > 2 * w:
                        # if there is a mutation close to hotspot
                        # despite LoF it is likely a driver and the second copy of this gene is broken
                        # this happens in sampl 905968  in PIK3CA gene (5290)
                        gene_scores[gene_id] = 1.0 / (dist + 1)
                    else:
                        # otherwise we set -1 to differentiate it from neutral mutations
                        # because oncogene might turn on after domain truncation
                        gene_scores[gene_id] = -1
                else:
                    # non-hotspot mutation in non-LoF oncogene - small positive value
                    gene_scores[gene_id] = max(1.0 / (dist + 1), w)
            else:
                # only LoF mutations, oncogene is likely not functional or
                # it is turned-on after truncting mutations
                # therefore we set -1 to differentiate it from neutral mutations
                gene_scores[gene_id] = -1


        else:  # if unknown_significance
            if "truncating_mutation" in mut_functions:
                gene_scores[gene_id] = -2.0 * w
            else:
                gene_scores[gene_id] = -w

    return pd.Series(gene_scores)


# ### Mutation / Hotspot Reference Sources
#
# Here we are trying to understand the volatility of the gene by the mutations that are within them. These are outlined in the data files from other data sources similar to OncoKB

# In[12]:
print(7)

hotspots2D_subst = pd.read_csv( data+"SNV_hotspots.csv")
hotspots2D_indels = pd.read_csv( data + "INDEL_hotspots.csv")
if verbose:
    print("1D mutations", hotspots2D_subst.shape[0], hotspots2D_indels.shape[0])
hotspots2D_subst = hotspots2D_subst.loc[~hotspots2D_subst["Amino_Acid_Position"].str.contains("splice")]
hotspots2D_subst = hotspots2D_subst[["Hugo_Symbol", "Amino_Acid_Position"]]
hotspots2D_subst.drop_duplicates(inplace=True)
hotspots2D_indels = hotspots2D_indels[["Hugo_Symbol", "Amino_Acid_Position"]]
hotspots2D_indels.drop_duplicates(inplace=True)
if verbose:
    print("1D hotspots", hotspots2D_subst.shape[0], hotspots2D_indels.shape[0])
cBioPortal_muts = pd.concat([hotspots2D_subst, hotspots2D_indels], axis=0)
cBioPortal_muts_ = cBioPortal_muts.loc[cBioPortal_muts["Amino_Acid_Position"].str.contains("-"), :]
cBioPortal_muts = cBioPortal_muts.loc[~cBioPortal_muts["Amino_Acid_Position"].str.contains("-"), :]
cBioPortal_muts["Amino_Acid_Position"] = cBioPortal_muts["Amino_Acid_Position"].apply(int)
# cBioPortal_muts["start"] = cBioPortal_muts["Amino_Acid_Position"]
# cBioPortal_muts["end"] = cBioPortal_muts["Amino_Acid_Position"]
cBioPortal_muts["aa_pos"] = cBioPortal_muts["Amino_Acid_Position"].apply(lambda x: (x, x))

# cBioPortal_muts_ ["start"] = cBioPortal_muts_["Amino_Acid_Position"].apply(lambda x: int(x.split("-")[0]))
# cBioPortal_muts_["end"] = cBioPortal_muts_["Amino_Acid_Position"].apply(lambda x: int(x.split("-")[1]))
cBioPortal_muts_["aa_pos"] = cBioPortal_muts_["Amino_Acid_Position"].apply(lambda x: tuple(map(int, x.split("-"))))

cBioPortal_muts = pd.concat([cBioPortal_muts_, cBioPortal_muts], axis=0)
cBioPortal_muts["source"] = "1D"

### 3D

cBioPortal_3D = pd.read_csv( data + "S5_3d_hotspots.csv")
cBioPortal_3D = cBioPortal_3D[["Gene", "Amino_Acid_Position"]]
if verbose:
    print("3D mutations", cBioPortal_3D.shape[0])
cBioPortal_3D.drop_duplicates(inplace=True)
if verbose:
    print("3D hotspots ", cBioPortal_3D.shape[0])
cBioPortal_3D["Amino_Acid_Position"] = cBioPortal_3D["Amino_Acid_Position"].apply(int)
# cBioPortal_3D["start"] = cBioPortal_3D["Amino_Acid_Position"]
# cBioPortal_3D["end"] = cBioPortal_3D["Amino_Acid_Position"]
cBioPortal_3D["aa_pos"] = cBioPortal_3D["Amino_Acid_Position"].apply(lambda x: (x, x))
cBioPortal_3D.rename({"Gene": 'Hugo_Symbol'}, axis="columns", inplace=True)
cBioPortal_3D["source"] = "3D"
cBioPortal_muts = cBioPortal_muts[["Hugo_Symbol", "aa_pos"]]
cBioPortal_3D = cBioPortal_3D[["Hugo_Symbol", "aa_pos"]]
cBioPortal_muts = pd.concat([cBioPortal_3D, cBioPortal_muts], axis=0)
cBioPortal_muts.sort_values(["Hugo_Symbol"], ascending=True, inplace=True)

cBioPortal_genes = set(cBioPortal_muts["Hugo_Symbol"].values)
if verbose:
    print("genes with hotspot", len(cBioPortal_genes))
cBioPortal_muts = cBioPortal_muts[["Hugo_Symbol", "aa_pos"]]
cBioPortal_muts.drop_duplicates(inplace=True)
if verbose:
    print("unique hotspots", cBioPortal_muts.shape[0])
# cBioPortal_muts


# #### Compare mutaitons with onco genes

# In[8]:
print(8)

if verbose:
    print("genes with and without mutations")
for kind in ["OG", "TSG", "OG_TSG"]:
    oncokb = OncoKB[kind].intersection(cBioPortal_genes)
    b = bailye[kind].intersection(cBioPortal_genes)
    if verbose:
       print(kind, "oncokb", len(oncokb), "/", len(OncoKB[kind]),
          "Bailye_2018", len(b), "/", len(bailye[kind]))

# In[9]:
print(9)

for kind in ["OG", "TSG", "OG_TSG"]:
    oncokb = OncoKB[kind]  # .intersection(cBioPortal_genes)
    b = bailye[kind]  # .intersection(cBioPortal_genes)
    if verbose:
        print(kind, "oncokb", len(oncokb), \
          "Bailye_2018", len(b), "union", len(b.union(oncokb)), \
          "mutated:", \
          len(b.union(oncokb).intersection(cBioPortal_genes)))

# In[14]:
print(10)

# TSG - not necessarity with a hotspot, take union of OncoKB and  Bailye_2018
TSG = OncoKB["TSG"].union(bailye["TSG"])
# OG - must have a hotspot
OG = OncoKB["OG"].union(bailye["OG"]).intersection(cBioPortal_genes)
if verbose:
    print("TSG:", len(TSG), "with hotspots:", len(TSG.intersection(cBioPortal_genes)),
      "\nOG with hotspots:", len(OG))
# list of hotspots:
hotspots = cBioPortal_muts.loc[cBioPortal_muts["Hugo_Symbol"].isin(OG | TSG), :]
if verbose:
    print("Hotspots (all from cBioPortal):", cBioPortal_muts.shape[0],
      "\n\ton OG:", hotspots.loc[hotspots["Hugo_Symbol"].isin(OG)].shape[0],
      "\n\ton TSG:", hotspots.loc[hotspots["Hugo_Symbol"].isin(TSG)].shape[0])
# hotspots.head(5)


# ### Often gene come with multiple names - This attempts to clean this up
#
# At the end we will have all genes translated to have
#
# Entrez Gene Id
#
# ```Entrez Gene maintains records from genomes which have been completely sequenced, which have an active research community to submit gene-specific information, or which are scheduled for intense sequence analysis```
#
# Hugo Symbol
#
# ``` The HUGO Gene Nomenclature Committee is a committee of the Human Genome Organisation that sets the standards for human gene nomenclature. ```
#
#
#
#

# In[13]:
print(11)

ncbi_synonyms = expand(NCBI[["Synonyms", "GeneID"]], column="Synonyms", sep="|")
ncbi_synonyms = parse_mapping_table(ncbi_synonyms, "Synonyms", "GeneID")

# In[16]:
print(12)

tsg_og_gene_names = pd.DataFrame.from_dict({"genes":dict(zip(list(OG | TSG),list(OG | TSG)))})
tsg_og_gene_names, query2target, not_mapped = apply_mappers(tsg_og_gene_names, ncbi_symbols, ncbi_synonyms, verbose = True,handle_duplicates = "keep")
OG_entrez = set(map(lambda x : query2target[x],OG))
TSG_entrez = set(map(lambda x : query2target[x],TSG))
#print(query2target)
hotspots["Entrez_gene_ID"] = hotspots["Hugo_Symbol"].apply(lambda x : query2target[x])

#hotspots.head()

# # Section 2
# ## Feature Mapping
#
# All the code above preps out look up tables and reference information. Now, we need to run all the variants from the WES variants from the MOLI (GDSC dataset)
#
# We will need to run the parts below against each patient WES file
#

# ## Load entire WES feature set used by MOLI paper (GDSC dataset)

# In[19]:
print(13)

wes_df = pd.read_csv(data + 'WES_variants.csv')
#wes_df.head()

# In[20]:
print(14)

mapped_wes_df = wes_df[["COSMIC_ID", "Gene", "AA", "Classification"]]
#mapped_wes_df.head(3)

# In[22]:
print(15)

# define functional effect of a mutation
gdsc_truncating = ['frameshift', 'nonsense', 'stop_lost', 'ess_splice']
gdsc_point = ['missense', 'Missense', 'inframe']
mapped_wes_df.loc[mapped_wes_df["Classification"].isin(gdsc_truncating), "function"] = "truncating_mutation"
mapped_wes_df.loc[mapped_wes_df["Classification"].isin(gdsc_point), "function"] = "point_mutation"
mapped_wes_df = mapped_wes_df.loc[mapped_wes_df["Classification"].isin(gdsc_truncating + gdsc_point), :]

plt.figure(figsize=(20, 5))
plt.subplot(1, 2, 1)
n_mutations_gdsc = mapped_wes_df.groupby("COSMIC_ID").size().sort_values(ascending=False)
tmp = plt.hist(n_mutations_gdsc, bins=50, range=(0, 2000))
tmp = plt.title("number of mutations per sample in GDSC")
ax = plt.subplot(1, 2, 2)
n_LoF = mapped_wes_df.loc[mapped_wes_df["Classification"].isin(gdsc_truncating), :].groupby(
    "COSMIC_ID").size().sort_values(ascending=False)
n_GoF = mapped_wes_df.loc[mapped_wes_df["Classification"].isin(gdsc_point), :].groupby("COSMIC_ID").size().sort_values(
    ascending=False)
df = pd.DataFrame.from_dict({"truncating_mutations": n_LoF, "point_mutations": n_GoF})
# tmp = plt.hist(n_LoF, bins=50)
# tmp = plt.title("LoF mutations per sample in GDSC")
#df.plot.hist(stacked=True, bins=50, ax=ax, range=(0, 2000))
if verbose:
    print("mean", round(np.mean(n_mutations_gdsc), 1), "median", np.median(n_mutations_gdsc))


# ### Helper funcs

# In[23]:
print(16)

def binarize(x):
    if x == 0:
        return x
    else:
        return 1


def get_aa_pos_GDSC(x):
    l = re.findall(r'\d+', x)
    if len(l) == 0:
        # print(x, "positions not found.",file = sys.stderr)
        return ("NA", "NA")
    elif len(l) == 1:
        return (int(l[0]), int(l[0]))
    else:
        return (int(l[0]), int(l[1]))


# ### Get actual position of mutations

# In[24]:
print(17)

mapped_wes_df["aa_pos"] = mapped_wes_df["AA"].apply(lambda x: get_aa_pos_GDSC(x))
#mapped_wes_df.head(5)

# ### Look to see if mutaitons are in onco genes

# In[25]:
print(18)

ensg_genes = []
symbols = []
for gene in set(mapped_wes_df["Gene"].values):
    if gene.startswith("ENSG"):
        ensg_genes.append(gene)

    else:
        symbols.append(gene)
if verbose:
    print("ENSEMBL gene IDs found:", len(ensg_genes))
ensg_genes = mapped_wes_df.loc[mapped_wes_df["Gene"].isin(ensg_genes), :]
gdsc_hugo = mapped_wes_df.loc[mapped_wes_df["Gene"].isin(symbols), :]
if verbose:
    print("Gene symbols:", len(set(gdsc_hugo["Gene"].values)))

# ### Map to Entrez

# In[26]:
print(19)

# mapping to Entrez
df = gdsc_hugo[["Gene"]].copy()
df.set_index("Gene", drop=False, inplace=True)
df.drop_duplicates(inplace=True)
df, query2target, not_mapped = apply_mappers(df, ncbi_symbols, ncbi_synonyms, verbose=True, handle_duplicates="keep")

# ### Map out to get entrez_gene_id
#
# ### This value, entrez_gene_id is important, it will become the name of all of our feature columns
#

# In[27]:
print(20)
gdsc_hugo = gdsc_hugo.loc[~gdsc_hugo["Gene"].isin(not_mapped), :].copy()
gdsc_hugo["Entrez_gene_ID"] = gdsc_hugo["Gene"].apply(lambda x: query2target[x])
#gdsc_hugo.head(10)

# ### Score each variant / mutation

# In[28]:
print(21)

w = 0.05
#print("THIS IS 2",TSG_entrez, OG_entrez, hotspots)
gdsc_scores = {}
for sample in list(set(gdsc_hugo["COSMIC_ID"].values)):
    df = gdsc_hugo.loc[gdsc_hugo["COSMIC_ID"] == sample, :]
    gdsc_scores[sample] = mutations2score(df, TSG_entrez, OG_entrez, hotspots, w=w)
gdsc_scores = pd.DataFrame.from_dict(gdsc_scores)
gdsc_scores.fillna(0, inplace=True)
#print("THIS IS 3")
# ### Save file with contious value scoring

# In[29]:
print(22)

gdsc_scores.to_csv(data + "GDSC.non_bin_mutations.tsv", sep="\t")
gdsc_scores.loc[[7157, 5290, 3845, 4893], :]
#print("THIS IS 4")
# ### Make decsion if we want to binarize the score

# In[30]:
print(23)

gdsc_scores = gdsc_scores.applymap(lambda x: binarize(x))
gdsc_scores.to_csv(data + "GDSC.binary_mutations.tsv", sep="\t")
gdsc_scores.loc[[7157, 5290, 3845, 2312, 4893], :]

# ### Prep data for transpose

# In[33]:
print(24)

# load one back in
cosmic_wes_binary = pd.read_csv(data + "GDSC.binary_mutations.tsv", sep="\t")
cosmic_wes_binary.set_index("Unnamed: 0", drop=True, inplace=True)
cosmic_wes_binary.index.names = ['Gene_ID']
#cosmic_wes_binary.head()

# In[34]:
print(25)

cosmic_wes_binary_transposed = cosmic_wes_binary.transpose()
#cosmic_wes_binary_transposed.head()

# ## Save out file transposed with COSMIC ID as rows

# In[37]:
print(26)

cosmic_wes_binary_transposed.index.names = ['COSMIC_ID']
cosmic_wes_binary_transposed.to_csv(data + "wes_binary_scored_transposed.tsv", sep="\t")

# In[38]:
print("Completed Wes Processing")

#cosmic_wes_binary_transposed.head()


