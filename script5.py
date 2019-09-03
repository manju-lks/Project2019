#!/usr/bin/env python
# coding: utf-8

# In[1]:


import vcf
import numpy as np
import pandas as pd
from collections import Counter


# In[2]:


vcf_sample1_snp = vcf.Reader(open("sample1_final_filtered_snps.vcf","r"))
vcf_sample2_snp = vcf.Reader(open("sample2_final_filtered_snps.vcf","r"))
vcf_sample3_snp = vcf.Reader(open("sample3_final_filtered_snps.vcf","r"))
vcf_sample4_snp = vcf.Reader(open("sample4_final_filtered_snps.vcf","r"))
vcf_sample5_snp = vcf.Reader(open("sample5_final_filtered_snps.vcf","r"))
vcf_sample6_snp = vcf.Reader(open("sample6_final_filtered_snps.vcf","r"))
vcf_sample7_snp = vcf.Reader(open("sample7_final_filtered_snps.vcf","r"))


# In[3]:


lis_vcfs=[vcf_sample1_snp, vcf_sample2_snp, vcf_sample3_snp, vcf_sample4_snp, vcf_sample5_snp,
          vcf_sample6_snp, vcf_sample7_snp]


# In[4]:


n = len(lis_vcfs)
count=0
res_list=[]
df=pd.DataFrame()
for vcff in lis_vcfs:
    count+=1
    dictionary={}
    sample=[]
    for record in vcff:
        sample.append([record.CHROM, record.POS, record.REF,' '.join([str(elem) for elem in (record.ALT)])])
    dictionary[count] = sample
    for val in dictionary.values():
        df[count] = pd.Series(val)
        for elements in val:
            res_list.append(elements)


# In[5]:


df


# In[6]:


df[1].head()


# In[7]:


df[2].head()


# In[8]:


res_list[:10]


# In[9]:


n=len(lis_vcfs)       
counts=Counter(tuple(x) for x in res_list)
val=[]
for k in counts.keys():
    if counts[k] == n:
        val.append(k)
        
print("number of snps common in all ", n, " files is ", len(val), " and they are, \n", val)


# In[ ]:




