import os
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import fisher_exact
from scipy.io import mmread
        
cl2ct = {4:'LMPP',
5:'CLP',
1:'HSC_MPP',
2:'MEP',
3:'CMP_BMP',
8:'GMP',
6:'Pro-B',
9:'MDP',
10:'pDC',
7:'Pre-B',
21:'Naive CD4 T1',
22:'Naive CD4 T2',
20:'Mature NK2',
17:'Basophil',
16:'Plasma cell',
19:'Mature NK1',
18:'Immature NK',
15:'Memory B',
11:'cDC',
13:'Monocyte 2',
12:'Monocyte 1',
14:'Naive B',
28:'Naive CD8 T3',
29:'Central memory CD8 T',
27:'Naive CD8 T2',
24:'Memory CD4 T',
23:'Naive Treg',
26:'Naive CD8 T1',
25:'Treg',
30:'Effector memory CD8 T',
31:'Gamma delta T'}

if __name__ == '__main__':
    
    df_barcode = pd.read_csv('./GSE129785_scATAC-Hematopoiesis-All.cell_barcodes.txt.gz', sep = '\t')
    id2fea = pd.read_csv('./GSE129785_scATAC-Hematopoiesis-All.peaks.txt.gz').Feature.values
    id2fea = [f.split('_') for f in id2fea]

    cluster = np.array([cl2ct[int(c[7:])] for c in df_barcode.Clusters.values])
    
    data_mat = mmread('./GSE129785_scATAC-Hematopoiesis-All.mtx')
    data_mat = data_mat.tocsr()
    
    for ct in set(cluster):
        ind = cluster == ct
        vec = np.array(data_mat[:, ind].mean(axis = 1))[:, 0]
        ind_peak = vec > 0.25
        diffpeaks = [id2fea[i] for i in np.where(ind_peak)[0]]
        print(ct, sum(ind_peak))
        pd.DataFrame(diffpeaks).to_csv('%s.bed' % ct.replace(' ', '_'), index = False, sep = '\t', header = False)