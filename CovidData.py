import numpy as np
import csv
import random
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
import seaborn as sns
import pandas as pd
df= pd.read_csv("/mnt/c/Biyoinformatik/Data Analysis/Data/GSE152418_p20047_Study1_RawCounts.txt", sep="\t", index_col=0) #datamız txt dosyasında olduğu için sütunlar arasında tab ile ayrılmıştır. index_col=0 ile ilk sütunu index olarak alıyoruz, sep="\t" ile de tab ile ayrıldığını belirtiyoruz.
#print(df.shape) #veri setimizin boyutunu yazdırıyoruz.
#print(df.head()) #veri setimizin ilk 5 satırını yazdırıyoruz.

df_filtered = df[df.sum(axis=1) >= 10] #veri setimizdeki toplam değeri 10'dan büyük olan satırları filtreliyoruz.
#print(df_filtered.shape) #filtrelenmiş veri setimizin boyutunu yazdırıyoruz.   
#axis=0 -> sütunlar üzerinde işlem yapar, axis=1 -> satırlar üzerinde işlem yapar. sum() fonksiyonu ile satırların toplamını alıyoruz ve >= 10 ile 10'dan büyük olanları filtreliyoruz.
#10 ve altındaki değerler, genellikle düşük ifade seviyelerine sahip genlerdir ve analizde gürültü oluşturabilir. Bu nedenle, bu tür genleri filtreleyerek daha güvenilir sonuçlar elde etmeyi amaçlıyoruz.
#CPM = (Counts Per Million) ifadesi, genlerin ifade seviyelerini normalize etmek için kullanılan bir yöntemdir. CPM, her genin ifade seviyesini toplam ifade seviyesine bölerek ve ardından 1 milyon ile çarparak hesaplanır. Bu sayede, farklı örnekler arasındaki ifade seviyeleri karşılaştırılabilir hale gelir.
#df_filtered -> her hücrede bir genin bir örnekteki ham sayımı
#df_filtered.sum(axis=0) -> her örneğin (sütunun) toplam sayımı (yani o örnekte kaç tane okuma olduğu, tüm genler toplanarak)

df_cpm = (df_filtered / df_filtered.sum(axis=0)) * 1e6 #veri setimizi CPM (Counts Per Million) formatına çeviriyoruz.
#df_filtered / df_filtered.sum(axis=0) -> her hücredeki değeri, o sütunun toplam değeri ile bölüyoruz. Bu sayede her hücredeki değer, o örnekteki toplam ifade seviyesine göre normalize edilmiş oluyor.
#print(df_cpm.shape) #CPM formatına çevrilmiş veri setimizin boyutunu yazdırıyoruz.
#print(df_cpm.head()) #CPM formatına çevrilmiş veri setimizin ilk 5 satırını yazdırıyoruz.
condition_labels = np.array([
    "Convalescent",  # S145_nCOV001_C
    "COVID-19",       # S147_nCoV001EUHM-Draw-1
    "COVID-19",       # S149_nCoV002EUHM-Draw-2
    "COVID-19",       # S150_nCoV003EUHM-Draw-1
    "COVID-19",       # S151_nCoV004EUHM-Draw-1
    "COVID-19",       # S152_nCoV006EUHM-Draw-1
    "COVID-19",       # S153_nCoV007EUHM-Draw-1
    "COVID-19",       # S154_nCoV0010EUHM-Draw-1
    "COVID-19",       # S155_nCOV021EUHM
    "COVID-19",       # S156_nCOV024EUHM-Draw-1
    "COVID-19",       # S157_nCOV0029EUHM
    "COVID-19",       # S175_nCoV024EUHM-Draw-2
    "COVID-19",       # S176_nCoV025EUHM-Draw-1
    "COVID-19",       # S177_nCoV025EUHM-Draw-2
    "COVID-19",       # S178_nCoV028EUHM-Draw-1
    "COVID-19",       # S179_nCoV033EUHM-Draw-1
    "COVID-19",       # S180_nCoV034EUHM-Draw-1
    "Healthy",        # S061_257
    "Healthy",        # S062_258
    "Healthy",        # S063_259
    "Healthy",        # S064_260
    "Healthy",        # S065_261
    "Healthy",        # S066_265
    "Healthy",        # S067_270
    "Healthy",        # S068_272
    "Healthy",        # S069_273
    "Healthy",        # S070_279
    "Healthy",        # S071_280
    "Healthy",        # S181_255
    "Healthy",        # S182_SHXA10
    "Healthy",        # S183_263
    "Healthy",        # S184_SHXA18
    "Healthy",        # S185_266
    "Healthy"         # S186_SHXA14
])
#print(list(df_filtered.columns)) #veri setimizin sütun isimlerini yazdırıyoruz.
mask = condition_labels != "Convalescent" #condition_labels değişkeninde "Convalescent" olan satırları filtreliyoruz.
df_cpm_filtered = df_cpm.loc[:, mask] #mask değişkeni ile filtrelenmiş sütunları alıyoruz.
condition_labels_filtered = condition_labels[mask] #mask değişkeni ile filtrelenmiş condition_labels değişkenini alıyoruz.

#print(mask.sum()) #kaç tane "Convalescent" olmayan örnek olduğunu yazdırıyoruz.
#print(df_cpm_filtered.shape) #filtrelenmiş veri setimizin boyutunu yazdırıyoruz.
#print(condition_labels_filtered) #filtrelenmiş condition_labels değişkenini yazdırıyoruz.

mean_control = df_cpm_filtered.loc[:, condition_labels_filtered == "Healthy"].mean(axis=1) #Healthy olan örneklerin ortalamasını alıyoruz.
mean_covid = df_cpm_filtered.loc[:, condition_labels_filtered == "COVID-19"].mean(axis=1) #COVID-19 olan örneklerin ortalamasını alıyoruz.
keep_express = (mean_control >= 5) | (mean_covid >= 5) #mean_control ve mean_covid değişkenlerinde 5'ten büyük olanları filtreliyoruz.
df_cpm_final = df_cpm_filtered.loc[keep_express, :] #keep_express değişkeni ile filtrelenmiş veri setimizi alıyoruz.
#print(df_cpm_final.shape) #filtrelenmiş veri setimizin boyutunu yazdırıyoruz.
#print(len(df_cpm_final.index)) #filtrelenmiş veri setimizin satır sayısını yazdırıyoruz.

from scipy.stats import ks_2samp

epsilon = 1 
p_values = []
log2_FCs = []

for i in range(df_cpm_final.shape[0]):
    gene = df_cpm_final.iloc[i, :]
    healty = gene.loc[condition_labels_filtered == "Healthy"]
    covid = gene.loc[condition_labels_filtered == "COVID-19"]
    ks_stat, p_val = ks_2samp(healty, covid)
    p_values.append(p_val)

    healthy_mean = healty.mean()
    covid_mean = covid.mean()
    FC = (covid_mean + epsilon) / (healthy_mean + epsilon)
    log2_FC = np.log2(FC)
    log2_FCs.append(log2_FC)

p_values = np.array(p_values)
log2_FCs = np.array(log2_FCs)

p_values_bonf = p_values * len(df_cpm_final)
deg_mask = (p_values_bonf < 0.05) & (np.abs(log2_FCs) > 1)
sig_genes = df_cpm_final.index[deg_mask]

sig_log2_FCs = log2_FCs[deg_mask]
sig_p_values_bonf = p_values_bonf[deg_mask]

sig_data = df_cpm_final.loc[deg_mask]
#print(f"Number of significant genes: {len(sig_genes)}")
#print(f"Total Number of Genes: {len(df_cpm_final)}")
#print(f"DEG percentage: {len(sig_genes) / len(df_cpm_final) * 100}")
#print()
#print(f"Number of upregulated genes: {np.sum(sig_log2_FCs > 0)}")
#print(f"Number of downregulated genes: {np.sum(sig_log2_FCs < 0)}")

#upregulated_genes = sig_genes[sig_log2_FCs > 0]
#for gene in upregulated_genes:
#       print(gene)

sig_df = pd.DataFrame({
    "Gene": sig_genes,
    "p_value_bonf": sig_p_values_bonf,
    "log2_FC": sig_log2_FCs
})

sig_data_log2fc = pd.DataFrame(
    np.zeros_like(sig_data),
    index=sig_data.index,
    columns=sig_data.columns
)
#print(sig_data_log2fc.shape) #sig_data_log2fc veri setimizin boyutunu yazdırıyoruz.

epsilon = 1
for i in range(len(sig_genes)):
    gene = sig_data.iloc[i, :]
    healthy_mean = gene.loc[condition_labels_filtered == "Healthy"].mean()
    sig_data_log2fc.iloc[i, :] = np.log2((gene + epsilon) / (healthy_mean + epsilon))

#sns_plot = sns.clustermap(
 #   sig_data_log2fc,
  #  xticklabels = sig_data.columns,
   # yticklabels = sig_genes,
    #cmap = "coolwarm",
    #vmin = -6, vmax = 6
#)
#sns_plot.ax_heatmap.set_xticklabels(sns_plot.ax_heatmap.get_xmajorticklabels(), fontsize = 2)
#sns_plot.ax_heatmap.set_yticklabels(sns_plot.ax_heatmap.get_ymajorticklabels(), fontsize = 4)

#sns_plot.savefig("CovidData_heatmap.png")
#plt.show()

upregulated_genes = sig_genes[sig_log2_FCs > 0]
downregulated_genes = sig_genes[sig_log2_FCs < 0]

#print("\n".join(upregulated_genes))
#print("---")
#print("\n".join(downregulated_genes))

import mygene

mg = mygene.MyGeneInfo()
result = mg.querymany(list(upregulated_genes), scopes='ensembl.gene', fields='symbol', species='human')

id_to_symbol = {r['query']: r.get('symbol', r['query']) for r in result}

upregulated_symbols = [id_to_symbol[gene_id] for gene_id in upregulated_genes]
#print("\n".join(upregulated_symbols))

# upregulated_symbols listesinde CDK1'in index'ini bul
idx = upregulated_symbols.index("CDK1")
gene_to_plot = upregulated_genes[idx]  # bu, CDK1'in Ensembl ID'si olacak
#print(gene_to_plot)

gene_to_plot = "ENSG00000170312"  # CDK1'in Ensembl ID'si

healthy_expressions = df_cpm_final.loc[gene_to_plot, condition_labels_filtered == "Healthy"]
covid_expressions = df_cpm_final.loc[gene_to_plot, condition_labels_filtered == "COVID-19"]

healthy_mean = healthy_expressions.mean()
healthy_sem = healthy_expressions.std() / np.sqrt(len(healthy_expressions))
covid_mean = covid_expressions.mean()
covid_sem = covid_expressions.std() / np.sqrt(len(covid_expressions))

conditions = ['Healthy', 'COVID-19']
means = [healthy_mean, covid_mean]
errors = [healthy_sem, covid_sem]

plt.figure()
#plt.bar(conditions, means, yerr=errors, capsize=5, color=['blue', 'red'])
#plt.ylabel('Expression Level (CPM)')
#plt.title('CDK1 Expression')
#plt.savefig("CDK1_expression_barplot.png")
#plt.show()

plt.figure()

def add_jitter(values, jitter_amount=0.05):
    return values + np.random.uniform(-jitter_amount, jitter_amount, size=len(values))

healthy_x = add_jitter(np.ones(len(healthy_expressions)) *1)
covid_x = add_jitter(np.ones(len(covid_expressions)) *2)

plt.scatter(healthy_x, healthy_expressions, color='blue', label='Healthy', s=10)
plt.scatter(covid_x, covid_expressions, color='red', label='COVID-19', s=10)

plt.errorbar(1, healthy_mean, yerr=healthy_sem, fmt='o', color='black', capsize=5)
plt.errorbar(2, covid_mean, yerr=covid_sem, fmt='o', color='black', capsize=5)

plt.xticks([1, 2], ['Healthy', 'COVID-19'])
plt.ylabel('Expression Level(CPM)')
plt.title('CDK1 Expression')
plt.legend()
plt.savefig("CDK1_expression_scatter.png")
plt.show()

plt.figure()
correlation_matrix = np.corrcoef(sig_data)
sns_plot = sns.clustermap(correlation_matrix,
                          xticklabels=sig_genes,
                          yticklabels=sig_genes,
                          cmap="coolwarm",
                          annot=False,
                          linewidths=.5,
                          figsize=(12, 10))

sns_plot.ax_heatmap.set_xticklabels(sns_plot.ax_heatmap.get_xmajorticklabels(), fontsize=4)
sns_plot.ax_heatmap.set_yticklabels(sns_plot.ax_heatmap.get_ymajorticklabels(), fontsize=4)

sns_plot.savefig("covid_degs_correlation_heatmap.png")
plt.show()
