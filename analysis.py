
import numpy as np
import csv
import random
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
import seaborn as sns
import pandas as pd
df = pd.read_csv("../Data/GSE150910_gene-level_count_file.csv", index_col=0)
#print(df.head())
samples = df.columns.tolist()  #Üst başlıklar -> Örnek İsimleri
data = df.values #İçerideki Tablo -> Saf sayısal Gen Sayımları
genes = df.index.tolist() # Sol Etiketler -> Gen İsimleri

#print("Number of Genes (Rows):", len(df.index))
#print("Number of samples (Columns):", len(df.columns))
#print("Full DataFrame Shape (Rows, Columns):", df.shape)
#pandas'ta .iloc ile tablonun yapısını bozmadan INteger Location kullanarak kesesir.
#print(df.iloc[0:5, 0:4])  # İlk 5 satır ve ilk 4 sütunu yazdırır
condition_labels = df.columns.str.split('_').str[0].values #Örnek isimlerinden koşul etiketlerini çıkarır ('chp_26' -> 'chp')
#.str.split('_') ile her bir örnek ismini '_' karakterine göre böler ve .str[0] ile ilk parçayı alır. Bu, koşul etiketlerini temsil eder.
#print(condition_labels[:5])  # İlk 5 koşul etiketini yazdırır
#print(np.unique(condition_labels))  # Benzersiz koşul etiketlerini yazdırır
#control_data = df.loc[:, condition_labels == 'chp']  # Kontrol grubu verilerini seçer
#print(pd.Series(condition_labels).value_counts())  
# Koşul etiketlerinin sayısını yazdırır. pd.Series pandas'ta tek sütunlu veri yapısını temsil eder ve .value_counts() ile her bir benzersiz değerin sayısını hesaplar.


mask = condition_labels != "chp"  # Kontrol grubu dışındaki örnekleri seçmek için bir maske oluşturur
df_filtered = df.loc[:, mask]  # Kontrol grubu dışındaki verileri seçer. mask eklediğimiz için koşulu yazmaya gerek kalmaz, direkt mak yazarız.
#df_filtered yaptıktan sonra kontrol amaçlı boyut doğrulama yapabiliriz.
condition_labels = condition_labels[mask]  # Filtrelenmiş koşul etiketlerini günceller
#print("Filtered DataFrame Shape):", df_filtered.shape)  # Filtrelenmiş veri çerçevesinin boyutunu yazdırır.
#print("Filtered Column Names Shape:", df_filtered.columns.shape)  # Filtrelenmiş sütun isimlerinin boyutunu yazdırır.
#print("Filtered Labels Shape:", condition_labels.shape)  # Filtrelenmiş etiketlerin boyutunu yazdırır.

df_cpm = df_filtered.div(df_filtered.sum(axis=0), axis=1) * 1e6  
# Filtrelenmiş verileri CPM (Counts Per Million) ile normalize eder. Her bir örneğin toplam sayısını alır, 1 milyon ile çarpar ve her bir genin sayısını bu toplamla böler.
#CPM yaparken aslında hastalar arasındaki okuma sayısı(sekans derinliği) farkını adilleştirmek/eşitlemek için yapıyoruz.
#Pandas'ta .sum() sayıları toplamaya yarar. Ama tablomuz 2 boyutlu (satır ve sütun) olduğu için axis parametresi ile hangi eksende toplama yapılacağını belirtmemiz gerekir. axis=0 satırları toplar, axis=1 sütunları toplar.
#.div() ise division yani bölme işlemi demektir. df.div(x) yazdığımızda df tablosundaki her bir elemanı x ile böler. axis=1 yazdığımızda ise her bir sütunun kendi toplamına bölünmesini sağlar.
#.sum() payda kısmında, .div() ise pay kısmında yer alır. Yani her bir genin sayısı, o örneğin toplam sayısına bölünür ve 1 milyon ile çarpılır.
#Eğer gen okumasını, o hastanın toplam okuma sayısına bölmezsek; Hasta B'deki bütün genler sırf cihaz o gün çok çalıştı diye "yüksek ifade edilmiş" gibi görünür.
#İşte bu yüzden her geni kndi hastasının toplamına bölüyoruz ki cihazın yarattığı adaletsizliği silip sadece genin gerçek oranına bakalım.

mean_control = df_cpm.loc[:, condition_labels == "control"].mean(axis=1)
mean_ipf = df_cpm.loc[:, condition_labels == "ipf"].mean(axis=1)
to_keep = (mean_control >= 5) | (mean_ipf >= 5)
df_cpm = df_cpm.loc[to_keep, :]
# .loc[:, condition_labels == "control"] virgülden önceki : sembolü tüm genleri (satırları) al demek
# virgülden sonraki condition_labels == "control" ise yalnızca kontrol grubuna ait olan 206 hastalık sütun süzgecini uygular.
# Pandas'ta axis=1 işlemi soldan sağa doğru satır boyunca yapılacağını söyler.
# Her bir genin (satırın) kontrol hastaları arasındaki ortalama CPM değerini tek bir sütun halinde hesaplar.
# to_keep = (mean_control >= 5) | (mean_ipf >= 5) işlemini yaparken biyolojik olarak bir genin anlamlı sayılması için en az bir grupta belirgin ifade edlimesi gerekir.
# | (veya) operatörü sayesinde, gen iki gruptan birinde 5 ortalamyı geçerse "TRUE" değerini alır.
#df_cpm.loc[to_keep, :] ile bu kez virgülden önce to_keep maskesini koyarak satırları süzeriz. Şartı sağlayan genler kalır ve virgülden sonraki : ile tüm sütunlar korunur.
#print(df_cpm.shape)  # Filtrelenmiş ve CPM normalize edilmiş veri çerçevesinin boyutunu yazdırır.
#print(len(df_cpm.index))  # Filtrelenmiş ve CPM normalize edilmiş gen sayısını yazdırır.

# etiket seçimi -> .loc
# indeksleme -> .iloc

first_gene = df_cpm.iloc[0, :]  # İlk genin tüm örneklerdeki CPM değerlerini alır
# .iloc[0, :] ile ilk satırdaki tüm sütunları seçeriz
control_gene_0 = first_gene.loc[condition_labels == "control"]  # İlk genin kontrol grubundaki CPM değerlerini alır
ipf_gene_0 = first_gene.loc[condition_labels == "ipf"]  # İlk genin ipf grubundaki CPM değerlerini alır
# df_cpm.iloc[0, :] doğrudan sayısal indeksle seçim yapar.
# 0 yazarak en üstteki ilk geni (0.satırı) alıyoruz. virgülden sonraki : ise tüm hastaları kapsamasını sağlar.
# first_gene.loc[condition_labels == "ipf"] ile aynı ilk gen için sadece ipf hastalarına ait cpm değerlerini süzer.
#print(control_gene_0.shape)
#print(ipf_gene_0.shape)

from scipy.stats import ks_2samp

ks_statistic, p_value = ks_2samp(control_gene_0, ipf_gene_0)  # İlk genin kontrol ve ipf grupları arasındaki Kolmogorov-Smirnov testini uygular
#ks_statistic(istatistik Değeri) iki grubun dağılım eğrileri arasındaki en büyük farkın mesafesidir. 0 ile 1 arasında değer alır; 1e ne kadar yakınsa iki grup birbirinden o kadar farklıdır.
#p_value: bu farkın şans eseri ortaya çıkma olasılığıdır. genelde 0.05'ten (veya daha küçük) düşük bir p-değeri, genin iki grup arasında istatikssel olarak anlamlı şekilde farklılaştığını söyler.

#Kolmogorov-Smirnov(KS) Testi, iki farklı grubun verilerinin aynı dağılımdan gelip gelmediğini kontrol eden görsel ve matematiksel bir testtir.
# İki grubun verilerini küçükten büyüğe dizer ve birikimli dağılım grafiklerini(CDF) üst üste koyar.
# İki grafik arasındaki en geniş dikay mesafeyi ölçer. Bu mesafeye ks_statistic(KS istatistiği) denir.
# Mesafe büyükse iki grup birbirinden belirgin şekilde farklılaşmıştır. Bu farkın tesadüf olup olmadığı ise p_value söyler. p<0.05 ise "Bu gen kontrol ve hasta grubunda istatistiksel olarak farklı davranıyor" deriz.

#Log2 Dönüşümü
#RNA-Seq verilerinde gen ifadeleri çok dengesizdir. Bir gen 2 okunurken diğeri 50k okunabilir.
# Ayrıca gen ifadeleri katlanarak (2 katı, 4 katı, 8 katı) büyür.
# Bu dönüşümü kullanmamızın nedeni aşırı uçları baskılayarak 50k gibi devasa sayıların istatistiği bozmasını engeller.
# Görselleştirmeyi kolaylaştırarak katlanışları eşit aralıklara dönüştürür.

#Epsilon nedir ve neden eklenir?
# Matematikte log2=0, yani tanımsızdır.
# Verimizde hiç ifade edilmemiş genler olabileceğinden, logaritma almadan önce bütün değerlere çok küçük sabit bit sayı(epsilon) eklenir.
# Biyoinformatikte genelde ortalığı bozmayacak kadar küçük bir değer örneğin $\epsilon = 1, veya $\epsilon = 10^-6) seçilir.
# Dönüştürülmüş Veri = log2(CPM+epsilon)

# En baştaki ilk geni (indeks 0) seçiyoruz.
first_gene = df_cpm.iloc [0, :]
#Bu genin değerlerini gruplara ayırıyoruz.
control = first_gene.loc[condition_labels == "control"]
ipf = first_gene.loc[condition_labels == "ipf"]
# df_cpm.iloc[0, :] -> df_cpm tablosunun ilk satırını (ilk geni) temsil eden tek boyutlu bir Pandas Serisi oluşturur.
# first_gene.loc[condition_labels == "control"] -> İlk genin 206 hasta arasındaki değerlerinde sadece kontrol hastalarına ait 103 tanesini "control" değişkenine atar.
# first_gene.loc[condition_labels == "ipf"] -> Aynı genin hasta grubuna ait 103 değerini "ipf" değişkenine atar.
# " == " -> eşitlik karşılaştırması (comparison) opertatörüdür. iki tarfınan birbirine eşit olup olmadığını sorgular ve geriye sadece "T" veya "F" döndürür.
ks_statistic, p_value = ks_2samp(control, ipf)
#print(f"ks_statistic: {ks_statistic}")
#print(f"p_value: {p_value}")
#ks_2samp(control, ipf) -> iki gruptaki değerleri karşılaştırır ve geriye 2 sayı döndürür.
#ks_statistic, p_value = ... -> ilk sayı ks_statistic, ikinci sayıyı p_value değişkenine atar.
# f"..." -> f-string, yani süslü parantez içindeki deperlerin değerini metin içine doğrudan koyup yazdırmayı sağlar.

#Sadece p-deperi bize iki grubun farklı olduğunu söyler ama genin ne kadar arttığını veya azaldığını söylemez. işte bunu ölçmek için FOLD CHANGE kullanıyoruz.
# Fold Change (Kat Değişimi) Nedir?
# Fold Change, hasta grubundaki gen ifadesi ortalmasının kontrol grubuna oranıdır:
# Fold Change= Expression Mean (IPF)/Expression Mean(Control)
# Eğer IPF ortalaması 10, COntrol ortalaması 5 ise -> Fold Change = 10/5=2 (Gen hastalarda 2 katına çıkmış).
# Eğer IPF ortalması 5, Control ortalaması 10 ise -> Fold Change = 5/10 = 0.5 (Gen hastalarda yarıya düşmüş).
# Log2 almamızın asıl nedeni simetri elde etmek. Örnek verdiğimiz eşitsizlikte 2 kat artış 2 ile, 2 kat azalış 0.5 ile gösteriliyor.
# log2 tabanına geçerek aslında bu düzensiz sayıları
# 2 kat artış -> log2(2)=+1
# Değişim yok -> log2(1)=0
# 2 kat azalış -> log2(0.5)=-1
# böylece 0'ın sağındaki artışlar +1, +2, +3, solundaki azalışlar ise simetrik olarak negatif -1, -2, -3 olur. 

#Kontrol grubunun ortalaması 0 ise, payda 0 olamayacağı için bilgisayar sıfıra bölme hatası (ZeroDivisionError) verir. Ayrıca log2(0) tanımsızdır.
#Bu yüzden ortalamaları alırken paya ve paydaya çok küçük bir sabit (örn;1) eklenir.
#log2(IPF Ortlaması + epsilon/Control Ortalaması+Epsilon)

epsilon = 1
#ortalamaları alıyoruz
control_mean = control.mean()
ipf_mean = ipf.mean()

#fold change ve log2 dönüşümü
FC = (ipf_mean + epsilon) / (control_mean + epsilon)
log2_FC = np.log2(FC)
#print(log2_FC)

#sonuçları toplayacağımız listeler
p_values = []
log2_FCs = []
epsilon = 1

#df_cpm DraFrame'i üzerinde döngü
for i in range(df_cpm.shape[0]): # range modülü 0'dan 11.*96'ya kadar sayı üretir, .shape ise bize df'deki toplam satır sayısını verir.
    #Bu döngüyle bilgisyara "0. satırdan başla, 11.096. satıra kadar her bir geni sırayla eline al ve aşağıdaki adımları tekrar et" diyoruz.
    
    # i. satırdaki geni Series olrak çekiyoruz.
    gene = df_cpm.iloc[i, :] # .iloc[i, :] komutu, i. indeksteki satırın (yani o anki genin) tüm sütunlardaki değerlerini pandas serisi olarak çeker.
    
    #pandas loc ile gruplara ayırıyoruz
    control = gene.loc[condition_labels == "control"] # .loc kullanarak etiketi "control" olan sütunları cntrol değişkenine
    ipf = gene.loc[condition_labels == "ipf"] # "ipf" olanları ise ipf değişkenine ayırıyoruz. BÖylece o genin sağlıklı ve hasta gruptaki ifade değerleini ikiye bölmüş oluyoruz.

    #KS testi p-değeri
    ks_stat, p_val = ks_2samp(control, ipf) # iki grubun dağılımını KS testine sokarç Fonksiyon iki değer döndürür, biz sadece p-değerini almak için "p-val" deşkenini tutarız.
    p_values.append(p_val)

    #Ortalamalar ve log2 Fold Change (Pandas .mean())
    c_mean = control.mean() #grupların ortalamalarını alır
    i_mean = ipf.mean()

    FC = (i_mean + epsilon) / (c_mean + epsilon) # sıfıra bölünme riski ortadan kaldırır.
    log2_fc = np.log2(FC)
    log2_FCs.append(log2_fc) #append komutu, hesaplanan p-değerini döngünün baında açtığımız boş p_values listesinin en sonuna bir eleman olarak ekler(böylece 1. genin sonucu 1. sıraya, 2. genin sonucu 2. sıraya dizlir.)
# Sonuçları diziye dönüştürüyoruz.
p_values = np.array(p_values)
log2_FCs = np.array(log2_FCs)

#Bonferroni Düzeltmesi: Her p-değerini toplam gen sayısı ile çarpıyoruz
p_values_bonf = p_values * len(df_cpm)

#to_keep: hem istatiksel hem biyolojik anlamlılık şartını sağlayan genler
to_keep = (p_values_bonf <= 0.05) & (np.abs(log2_FCs) >=2)

#Gen isimlerini DataFrame index'inden çekiyoruz(sadece anlamlı genlerin isimleri)
sig_genes = df_cpm.index[to_keep]

#ilgili istatistikleri filtreliyoruz
sig_log2_FCs = log2_FCs[to_keep] #(anlamlı genlerin kat değişimi)
sig_p_values_bonf = p_values_bonf[to_keep] #sadece anlamlı genlerin p-değeri

#Orijinal CPM verisinin sadece anlamlı genlere ait satırlarını çekiyoruz
sig_data = df_cpm.loc[to_keep] #sadece anlamlı genlerin tüm CPM verisi

#print(f"Number of DEGs: {len(sig_genes)}") # filtreden geçen (hem istatiksel hem biyolojik olarak anlamlı gen sayısı)
#print(f"Total Number of Genes: {len(df_cpm)}") # toplam test ettiğim gen sayısı
#print(f"DEG Percentage: {len(sig_genes) / len(df_cpm) * 100}") #yüzde kaçının anlamlı çıktığı
#print()
#print(f"Number of Upregulated DEGs: {np.sum(sig_log2_FCs > 0)}") #boolean dizisi verrir, her eleman için "bu genin log2FV'si pozitif mi?" sorusunun cevabı.
#print(f"NUmber of Downregulated DEGs: {np.sum(sig_log2_FCs < 0)}") #negatif log2FV'ye sahip genler
#log2FC pozitifse, IPF grubundaki ortalama control grubundan büyük demektir. Yani gen hastalıkta daha fazla ifade edilmiş. Negatifse tam tersi, gen hastalıkta baskılanmış.
#Bu bilgiler, "hangi genler hastalıkta artmış, hangileri azalmış" sorusuna sayılsal bir cevap vermiş oluyoruz. Bu genelde rapotun ilk özet cümlesi olur ("79 DEG bulundu, 65'i uprefulated, 14'ü downregulated").

upregulated_genes = sig_genes[sig_log2_FCs > 0]
for gene in upregulated_genes:
    #print(gene)

#Bu hareketimizden sonra elimizde bir gen listesi oluşuyor. Bu gen listesini "Enrichr" sitesinde arattığımızda "bu genler hangi bioylojik yolaklar/hastalıklarla ilişkili" diye sorgulama aracı.
#Yani DEG listesi çıktıktan sonra genelde biyoinformatikte atılan bir sonraki adım, ham istatistikten biyolojik anlam çıkarmaya geçiş.
# Bu sitede araştırma yaptıktan sonra asıl bakmamız gerekenler "Pathways", "Ontologies(GO Biological Process)", "Diseases/Grugs" kısımlarına bakmamız gerekiyor.
# Pathways -> Burada ipf hastalığı kollajenle ilgili olduğu için "collagen biosynthesis" gibi terimler bu hastalıkla doğrudan örtüşür.
# Ontologies(GO Biological Process) -> genlerin hangi biyolojik süreçlere karıştığını gösterir.
# Disieases/Drugs -> genlerin hangi hastalıklarla ilişkilendirildiğini gösterir. Bu muhtemel en çarpıcı sonucu vereceği sekmedir. 


#import pandas as pd
#import os

# Önce tüm sonuçları tek bir DataFrame'de topluyoruz (isteğe bağlı ama pratik)
 sig_df = pd.DataFrame({
    "gene": sig_genes,
    "p_value": sig_p_values_bonf,
    "log2_FC": sig_log2_FCs
})

# Klasör yoksa oluştur
#os.makedirs("output", exist_ok=True)

# Upregulated (log2FC > 0) olanları filtreleyip kaydediyoruz
#sig_df[sig_df["log2_FC"] > 0].to_csv("output/upregulated_DEGs.csv", index=False)

# Downregulated (log2FC < 0) olanları filtreleyip kaydediyoruz
#sig_df[sig_df["log2_FC"] < 0].to_csv("output/downregulated_DEGs.csv", index=False)

#Bu adımda Heatmap hazırlıyoruz, bununla birlikte her hücrenin aynı ölçekte (log2FC bazında) renklenmesi için boş bir matris hazırlıyoruz.
sig_data_log2fc = pd.DataFrame(
    np.zeros_like(sig_data), #sig_data ie aynı boyutta (79genx206örnek), içi tamamen sıfırlarla dolu bir dizin oluşturur. Burada kalıp oluşturuyoruz.
    index = sig_data.index,
    columns = sig_data.columns
)
#print(sig_data.shape) #79, 206 demek aslında, anlamlı 79 genin (DEG), her biri için 206 örnekte (control+IPF hastalarının toplamı) CPM deperi olduğu anlamına gelir.

epsilon = 1 # sıfıra bölünmeyi önlemek için küçük bir sabit
for i in range(len(sig_genes)):
   gene = sig_data.iloc[i, :] #i. genin tüm örneklerdeki (206 örnek) CPM değerlerini bir pandas serisi olarak çekiyoruz.
   control_mean = gene.loc[condition_labels == "control"].mean() # bu genin sadece control grubundaki örnkelerinin ortalamasını alıyoruz -> bu, "referans noktası" oluyor.
   sig_data_log2fc.iloc[i, :] = np.log2((gene + epsilon) / (control_mean + epsilon)) #buradaki kritik fark şu: gene tek bir sayı değil, 206 elemanlı bir dizi/seri. yani bu satır, genin her bir örnekteki
# değreii tek tek control_mean'e bölüp log2 alıyor - control örnekleri de dahil olmak üzere (control örnekleri kendi ortalamasına bölündüğü için 0'a yakın çıkacak, IPF örnrkleri ise farklılık gösterecek).
# sig_data_log2fc.iloc[i, :] = ... -> hesaplanan bu 206 elemanlı diziyi, boş matrisimizin i. satırına yerleştiriyoruz.

#Heatmap'te her hücrenin (gen x örnek) rengi, o örneğin control ortalamasına göre ne kadar sapma gösterdiğini yansıtacak.
#Control örnkleri genelde 0'a yakın (yani "değişim yok") renkte görünecek, IPF örnekleri ise pozitif/negatif sapma gösterecek - bu da hastalık grubundaki değişimi görsel olarak net bir şekilde ortaya koyacak.

#Çalışmanın devamında aynı kodlara ihtiyacımız olacağı için heatmap kodlarını yoruma çevirdik.
#sns_plot = sns.clustermap(
 #  sig_data_log2fc,
  # xticklabels = sig_data.columns, #örnek isimleri (samples yerine)
  # yticklabels = sig_genes,        #gen isimleri
   #cmap = "coolwarm",
   #vmin = -6, vmax=6
#)
# x ve y eksenindeki yazı boyutlarını ayarlıyoruz.
#sns_plot.ax_heatmap.set_xticklabels(sns_plot.ax_heatmap.get_xmajorticklabels(), fontsize = 2)
#sns_plot.ax_heatmap.set_yticklabels(sns_plot.ax_heatmap.get_ymajorticklabels(), fontsize = 4)

#sns_plot.savefig("heatmap_for_analysis.pdf")
#plt.show()

# ======== GRAFİK OKUMA ========
# Renk Skalası (-6 ile +6 arası): Kırmızı -> o örnekte genin ifadesi control ortalamasına göre artmış(upregulated),
# mavi -> azalmış(downregulated), beyaz/açıkton -> control'e yakın(değişim yok).
#Üstteki dendrogram (ağaç yapısı, sütunlar için): Örnekleri (206 tanesi) birbirine ne kadar benzediklerine göre kümeliyor. GÖdrüğün gibi net iki büyük grup oluşuyor - bu muhtemelen control ve IPF örnkelrinin kendi aralarında kümlenmesi.
#Bu analizin doğru çalıştığının görsel bir kanıtı: EĞER HASTALIK GERÇEKTEN BU GENLERDE TUTARLI BİR ETKİ YARATIYORSA; ÖRNEKLER OTOMARİK OLARAK "HASTA" ve "SAĞLIKLI" DİYE İKİYE AYRILMALI
#Soldaki Dendrogram (satırlar, genler için): 79 geni, ifade patenrelerine göre kümeliyor. Üstteki büyük kırmızı blok, IPF grubunda birlikte artan genleri gösteriyor; alttaki mavi blok ise IPF'de birlikte azalan genleri gösteriyor.
#vmin= -6, vmax=6 neden önemli: Bu, renk skalasını tüm matriste sabitliyor- yani her hücre aynı cetvel ile ölçülüyor, bu da genler arası karşılaştırmayı adil hale getiriyor.
#Sonuç olarak bu görsel, tüm istatikseli analizinin (p-value, Bonferroni, log2fc eşikleri) gerçek bir biyolojik ayrımı (hasta vs sağlıklı) yakalanıldığını net bir şekilde gösteriyor.

#print(sig_data_log2fc.iloc[0, :10])
#print(sig_data_log2fc.iloc[5, :10])
#print(sig_data_log2fc.iloc[20, :10])

#Sırada tek bir genin control ve IPF grypları arasındaki ortlama ifade düzeyini, hata çubuklarıyla (SEM-Standart Error of Mean) birlikte gösteren bir bar grafiği

#gene_to_plot = "POSTN"

#control_expression = df_cpm.loc[gene_to_plot, condition_labels == "control"] #df_cpm'in index'i zaten gen isimleri, yani .loc["POSTN", ...] direkt o genin satırını seçip, ikinci parametre olan boolen maskeyle
#(condition_labels == "control") sadecce control örneklerini filtreliyor. 
#ipf_expression = df_cpm.loc[gene_to_plot, condition_labels == "ipf"]

#control_mean = control_expression.mean() #Control grbundaki tüm örenkerin POSTN geni için ölçülen CPM değerlerinin aritmetik ortalaması.
#control_sem = control_expression.std() / np.sqrt(len(control_expression)) #.std() standart sapma-standard deviation -> grup içindeki örnekleri birbirinden ne kadar farklı/dağınık olduğunu ölçer. Eğer tüm control örnekleri birbirine çok yakın değerler veriyorsa, std düşük olur.
#eğer örnekler araıdnaki büyük farklar varsa std yüksek olur.
#ipf_mean = ipf_expression.mean()
#ipf_sem = ipf_expression.std() / np.sqrt(len(ipf_expression)) #sem gruplarındaki bu np.sqrt(len(control_expression)) (örnek sayısının karkökü): kaç tane örnek olduğunun karekökü
# Std'yi bu sayıya bölmek (yani SEM hesaplamak): Şunu soruyoruz aslında -> "Ortalaamız(mean) ne kadar güvenilir?"
# Mantık şu- ne kadar çok örneğin varsa, ortalaman o kadar güvenilir olur (çünkü tesadüfi sapmalar birbirini dengeler).
# Sem, bu "ortalamını ne kadar sağlam olduğu" bilgisini tek bir sayıya indirger.
#conditions = ['Control', 'IPF']
#means = [control_mean, ipf_mean]
#errors = [control_sem, ipf_sem]

#Neden std yerine sem kullanılıyor?
#std -> bireysel örneklerin ne kadar iyi dağıldığını gösterir.
#sem -> ortalamanın kensinin ne kadar güvenilir olduğunu gösterir.
#plt.bar(conditions, means, yerr = errors, capsize = 5, color = ['blue', 'red'])
#plt.ylabel('Expression Level (CPM)')
#plt.title(f'{gene_to_plot} Expression')
#plt.savefig(f"heatmap_for_analysis_{gene_to_plot}_expression_bar.pdf")
#plt.show()
#Grafikte gördüğün hata çubukları (error bars) işte bu SEM değerini gösteriyor. Küçük hata çubuğu -> "bu ortalamaya güvenebiliriz, örnekler tutarlıydı".
#Büyük hata çubuğu -> "örnekler arasında çok fark vardı, ortalam daha az güvenilir."
#BİLİMSEL OLARAK NEDEN ÖNEMLİ?
#Grafikte control ve IPF çubukları arasında fark olduğunu görüyoruz (mesela 90 vs 400). Ama hata çubukları da küçükse, bu farkın tesadüfi değil, gerçek ve güvenilir bir fark olduğunu düşünmeni sağlar.
#Eğer hata çubukları çok büyük olsaydı ve iki grup birbirine "değebilecek" kadar yakın olsaydı, o zaman "bu fark gerçek mi yoksa şans mı?" diye şüphelenmemiz gerekirdi.
# Bu kavram aslında confidence interval (güven aralığı) ve t-test gibi konuların temelini oluşturur.

#Bu adıma geçtiğimizde ise Scatter Plot (dağılım grafiğini) bulmaya çalışacağız. Bir öncesindeki bar plot'un gelişmiş halidir.
#Bar grafiği sadec ortalamayı gösterirken, scatter plot her bir örneği tek tek nokta olarak gösteriyor. Bu çok daha bilgilendirici çünkü ortalamaının arkasındaki gerçek veri dağılımını görebilyoruz.(kaç örnek var, ne kdar yayılmışlar, aykırı değerler var mı gibi).

#gene_to_plot = "POSTN"

#control_expression = df_cpm.loc[gene_to_plot, condition_labels == "control"]
#ipf_expression = df_cpm.loc[gene_to_plot, condition_labels == "ipf"]

#control_mean = control_expression.mean()
#control_sem = control_expression.std() / np.sqrt(len(control_expression))
#ipf_mean = ipf_expression.mean()
#ipf_sem = ipf_expression.std() / np.sqrt(len(ipf_expression))

#def add_jitter(values, jitter_amount = 0.05):
 #   return values + np.random.uniform(-jitter_amount, jitter_amount, size = len(values))
#add_jitter fonksiyonu -> normalde tüm control noktalarını x=1 üzerine, tüm IPF noktalarını x=2 üzerine koyunca noktalar üst üste binip tek bir dikey çizgi gibi görünüp,
#ayırt edilemez. jitter(titreşim), her noktaya küük rastgele bir yatay kaymalar (-0.05 ile +0.05 arası) ekleyerek noktaları yatayda hafifçe dağıtıyor-böylece binlerce nokta
#üst üste binmeden, "bulut" gibi görünür hale geliyor. Bu, veri görselleştirmede "strip plot" veya "jitter plot" olarak bilinen standart bir tekniktir.

#control_x = add_jitter(np.ones(len(control_expression)) * 1)
#ipf_x = add_jitter(np.ones(len(ipf_expression)) * 2)

#plt.scatter(control_x, control_expression, color = 'blue', label = 'Control', s = 10)
#plt.scatter(ipf_x, ipf_expression, color = 'red', label = 'IPF', s = 10)
#plt.errorbar(1, control_mean, yerr = control_sem, fmt = 'o', color = 'black', capsize = 5)
#plt.errorbar(2, ipf_mean, yerr = ipf_sem, fmt = 'o', color = 'black', capsize = 5)

#plt.xticks([1, 2], ['Control', 'IPF'])
#plt.ylabel('Expression Level (CPM)')
#plt.title(f'{gene_to_plot} Expression')
#plt.legend()
#plt.savefig(f"{gene_to_plot}_expression_scatter.pdf")
#plt.show()

#Bar grafiğinda sadece "ortalama 90 vs 400" görüyordum. Bu grafikte ise control grubunun çoğu örneğinin gerçekten düşük(0-200 CPM) kümeleştiğini,
#ama IPF grubunun çok daha geniş bir aralığa (0'dan 1200'e kadar) yayıldığını görüyorum- yani hastalık sadece ortalamayı yükseltmiyor,
# aynı zamanda hastalar arası değişkenliği de artırıyor(bazı hastalarda çok yüksek, bazılarında ise control'e yakın). 
#Bu, ortalam tek başına söyleyemeceği önemli bir detay.

correlation_matrix = np.corrcoef(sig_data)

sns_plot = sns.clustermap(correlation_matrix,
                          xticklabels = sig_genes,
                          yticklabels = sig_genes,
                          cmap = "coolwarm",
                          annot = False,
                          linewidths = .5,
                          figsize = (12, 10))

sns_plot.ax_heatmap.set_xticklabels(sns_plot.ax_heatmap.get_xmajorticklabels(), fontsize = 4)
sns_plot.ax_heatmap.set_xticklabels(sns_plot.ax_heatmap.get_ymajorticklabels(), fontsize = 4)

sns_plot.savefig("degs_correlation_heatmap.pdf")
plt.show()

#Grafikte kırmızı -> iki gen birlikte artıp azalıyor(pozitif korelasyon), mavi-> ters yönde hareket ediyorlar(negatif korelasyon),
#sol üstteki kırmızı blok kendi aralarında güçlü korelasyon gösteriyor, büyük orta-sağ blok da (kolejen genleri) kendi aralarında kümeleniyor
# yani genler rastgele değil, birlikte hareket eden gruplar(modeller) halinde davranıyor, bu da biyolojik olarka anlamlı(aynı yolağın parçası ola genler birlikte regüle olur).
