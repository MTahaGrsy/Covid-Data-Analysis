import csv
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import pandas as pd

np.random.seed(42) #reproducibility (tekrarlanabilirlik) için
#statsmodels ne işe yarar?
#scipy.stats (RNA-Seq projesinde ks_2samp için kullanılan) daha çok "iki grup farklı mı?" gibi basit tesrker için uygundu.
#statsmodels ise regresyon modelleri (lineer regresyon, lojistik regresyon gibi) kurmak için kullanılıyor- GWAS'ta her SNP (genetik varyant) için "Bu varyant, hastalık durumunu ne kadar iyi açıklıyor" sorusuna regresyon ile cevap aranılacak.
#o yüzden bu kütüphaneyi kullanıyoruz.

#np.random.seed(42) neden önemli?
#Bu, rstgele sayı üreticinin "başlangıç noktasını" sabitliyor. Mock veri oluştururken rastgelelik kullanılacağı için (gerçek veri olmadığından simüle ediyoruz), bu satır olmadan her çalıştığında farklı rastgele veriler üretilir. 
#seed(42) sayesinde, rastgele sayılar üretilir- böylece sonuçların birebir karşılaştırılabilir olur.

#yazacağımız kod -> Tek bir SNP için sahte (mock) veri üretiyor.
N = 200 # toplam kişi sayısı (100 control + 100 "zombi" vakası)

p_control = 0.2 # control grubunda "riskli allel"in frekansı
p_case = 0.8 # vaka (zombi) grubunda "riskli allel"in frekansı

#Genotip nedir, neden 0/1/2 kodlanıyor?
#GWAS'ta her kişinin bir SNP (genetik varyant) için ik kopyası vardır (biri anneden, biri babadan).
#Genotip şu şekilde kodlanır:
#0 -> riskli alleli hiç taşımıyor (iki kopya da "normal")
#1 -> bir kopya riskli, bir kopya normal (heterozigot)
#2 -> iki kopya da riskli (homozigot)

snp_control = np.random.choice([0,1,2], size = int(N/2), p = [(1-p_control) ** 2, 2 * p_control * (1-p_control), p_control ** 2])
snp_case = np.random.choice([0,1,2], size = int(N/2), p = [(1-p_case) ** 2, 2 * p_case * (1-p_case), p_case ** 2])
#yukarıdaki satırda, Hardy-Weinberg dengesi denen bir genetik prensibi kullanıyor- popülasyonda bir allel frekansı biliniyorsa, genotip frekanslarını matematiksel olarak tahmin edebilirsin: 
# Genotip 0 olma olasılığı: (1-p)^2 (iki "normal" allel bir araya gelme olasılığı)
# Genotip 1 olma olasılığı: 2p(1-p) (riskli +bir normal, iki farklı sırada olabileceği için 2 ile çarpılır)
# Genotip 2 olma olasılığı: p^2 (iki riskli allel bir araya gelme olasılığı)

# COntrol grubunda riskli allel frekanası düşük (p_control = 0.2) olduğu için, control grubundaki insanların çoğu genotip 0 (riskli alleli yok) çıkacak.
# Vaka (zombi) grubunda ise riskli alleli frekansı yüksek(p_case = 0.8) olduğu için, çoğu kişi genotip 2 çıkacak - yani bu SNP ile "zombi olma" arasında güçlü bir ilişki olacak şekilde bilerek kurgulanmış bir senaryo.
snp = np.concatenate([snp_control, snp_case]) #iki grubu birleştiriyoruz (200 kişilik tke dizi)
labels = np.array([0] * 100 + [1] * 100) # 0 = control, 1 = vaka (zombi)

#Bu blok, "gerçekte SNP ile hastalık arasında güçlü bir ilişki olsaydı, veri nasıl görünürdü" sorusuna cevap veren yapay bir örnek oluşturuyor.