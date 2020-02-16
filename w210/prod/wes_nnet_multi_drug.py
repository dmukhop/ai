
# coding: utf-8

# In[2]:


from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split


# Install TensorFlow

import tensorflow as tf
import tensorflow_hub as hub


# ## Data Load
# 
# Loading random 5 drugs. Potentially we ould do a set with all the data, or drugs targeted for a particular cancer.

# In[3]:


wes_efgr = pd.read_csv("../data/wes_ic50_Erlotinib_binary_resp.tsv",sep = "\t", index_col="COSMIC_ID")
wes_pacl = pd.read_csv("../data/wes_ic50_Paclitaxel_binary_resp.tsv",sep = "\t", index_col="COSMIC_ID")
wes_suni = pd.read_csv("../data/wes_ic50_Sunitinib_binary_resp.tsv",sep = "\t", index_col="COSMIC_ID")
wes_sora = pd.read_csv("../data/wes_ic50_Sorafenib_binary_resp.tsv",sep = "\t", index_col="COSMIC_ID")
wes_rapa = pd.read_csv("../data/wes_ic50_Rapamycin_binary_resp.tsv",sep = "\t", index_col="COSMIC_ID")


# In[4]:


wes_multi_drug = pd.concat([wes_efgr, wes_pacl, wes_suni, wes_sora,wes_rapa], axis=0)


# In[5]:


wes_efgr.describe()


# In[6]:


#single drug output
wes_efgr.head()


# In[7]:


wes_multi_drug.head()


# In[8]:


len(wes_multi_drug)


# ## Prep Data
# Convert label to int, drop float version of response

# In[9]:


wes_multi_drug.loc[(wes_multi_drug.BINARY_RESPONSE == 'S'),'BINARY_RESPONSE'] = 1
wes_multi_drug.loc[(wes_multi_drug.BINARY_RESPONSE == 'R'),'BINARY_RESPONSE'] = 0

wes_multi_drug = wes_multi_drug.drop('LN_IC50', axis=1)


# In[10]:


# split data
train_data_set, dev_data_set, train_labels, dev_labels = train_test_split(wes_multi_drug, wes_multi_drug['BINARY_RESPONSE'], test_size=0.20, random_state=0)


# ### Translate to TF Data

# In[11]:


train_dataset = tf.data.Dataset.from_tensor_slices((train_data_set.values, train_labels.values))
train_dataset = train_dataset.shuffle(len(train_data_set)).batch(1)


# In[12]:


test_dataset = tf.data.Dataset.from_tensor_slices((dev_data_set.values, dev_labels.values))
test_dataset = test_dataset.shuffle(len(dev_data_set)).batch(1)


# ## Build Simple Model

# In[13]:


model = tf.keras.models.Sequential([
  tf.keras.layers.Dense(units=128, activation='relu'),
  tf.keras.layers.Dropout(0.2),
  tf.keras.layers.Dense(10, activation='softmax')
])

model.compile(optimizer='adam',
              loss='sparse_categorical_crossentropy',
              metrics=['accuracy'])


# In[15]:


model.fit(train_dataset, epochs=10, steps_per_epoch=150)


# ## Eval Model

# In[16]:


model.evaluate(test_dataset, verbose=2, steps=len(dev_data_set))


# In[17]:


model.summary()


# In[18]:


model.save('wes_multi_model')


# ## Load Saved Model As Layer

# In[19]:


hub_url = "wes_multi_model/"
embed = hub.KerasLayer(hub_url, input_shape=(18384,), dtype=tf.float64, trainable=False)


# In[20]:


wes_multi_base_model = tf.keras.Sequential([
    embed,
    tf.keras.layers.Dense(16, activation="relu"),
    tf.keras.layers.Dense(1, activation="sigmoid"),
])


# In[21]:


wes_multi_base_model.compile(optimizer='adam',
              loss='binary_crossentropy',
              metrics=['accuracy'])


# In[22]:


wes_multi_base_model.fit(train_dataset, epochs=10, steps_per_epoch=150)


# In[23]:


wes_multi_base_model.evaluate(test_dataset, verbose=2, steps=len(dev_data_set))


# In[24]:


# wes_multi_base_model.summary()

