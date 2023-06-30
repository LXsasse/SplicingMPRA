import numpy as np
import sys, os

inputfeat = ['MFE_main', 'MFE_ens', 'Freq_of_main', 'Ensemble_div', 'GC_content', 
        'alt_SE', 'low_hex', 'top_hex', 'high_motif', 'low_motif', 'mod_motif', 'Yof_high_motif', 'Yof_low_motif', 'fiveP_ss', 'BP']

seqname = 'barcode' 
outputfeat = 'ave_SE'

def isfloat(f):
    try:
        f = float(f)
        return f
    except:
        f = float(f.replace('-', 'e-').replace('+', 'e+'))
        return f


def check(line):
    try:
        a= np.array(line)[columns].astype(float)
        return a
    except:
        a = []
        for i, fi in enumerate(line):
            if (('-' in fi[1:]) or ('+' in fi[1:])) and 'e' not in fi:
                fi = fi[0] + fi[1:].replace('-', 'e-').replace('+', 'e+')
            a.append(float(fi))
        return np.array(a)


data = open('complete_master_both_sets_ilib_ave.tsv', 'r').readlines()
ids, feats, ef, seqs = [], [], [], []
for i, line in enumerate(data):
    if i == 0:
        colnames = line.strip().split('\t')
        columns = np.isin(line.strip().split('\t'), inputfeat)
        out = line.strip().split('\t').index(outputfeat)
        sc = line.strip().split('\t').index('sequence')
        print(len(inputfeat), np.sum(columns))
    else:
        line = line.strip().split('\t')
        if len(line[0]) > 10:
            ids.append('seq_'+line[0])
            feat = np.array(line)[columns]
            feats.append(check(feat))
            ef.append(line[out])
            seqs.append(line[sc])

ids = np.array(ids)
ef = np.array(ef)
feats = np.array(feats)
seqs = np.array(seqs)

sort = np.argsort(ids)
ids, ef, feats, seqs = ids[sort], ef[sort], feats[sort], seqs[sort]



np.savetxt('Pred_meanefficiency.txt', np.array([ids, ef]).astype(str).T, fmt = '%s')
np.savez_compressed('SeqFeat_input.npz', names = ids, features = np.array(feats, dtype = float), featuretypes = np.array(colnames)[columns])

if '--kmers' in sys.argv:
    def generatek(l):
        kl = np.array(list('ACGT'))
        km = np.copy(kl)
        for i in range(l-1):
            kc = []
            for kt in km:
                for c in kl:
                    kc.append(kt+c)
            km = np.array(kc)
        return np.sort(km)

    for k in range(3,8):
        print(k)
        kmers = generatek(k)
        print(kmers)
        kfeat = np.zeros((len(seqs), len(kmers)), dtype = np.float16)
        for s, seq in enumerate(seqs):
            kset = []
            for l in range(len(seq)-k+1):
                kset.append(seq[l:l+k])
            uk, un = np.unique(kset, return_counts = True)
            kfeat[s, np.isin(kmers, uk)] = un
        np.savez_compressed('SeqFeatkmer'+str(k)+'_input.npz', names = ids, features = np.append(kfeat, np.array(feats, dtype = np.float32),axis = 1), featuretypes = np.append(kmers, np.array(colnames)[columns]))
        print('SeqFeatkmer'+str(k)+'_input.npz')
        np.savez_compressed('Seqkmer'+str(k)+'_input.npz', names = ids, features = kfeat, featuretypes = kmers)
        print('Seqkmer'+str(k)+'_input.npz')

if '--barcodekmers' in sys.argv:
    def generatek(l):
        kl = np.array(list('ACGT'))
        km = np.copy(kl)
        for i in range(l-1):
            kc = []
            for kt in km:
                for c in kl:
                    kc.append(kt+c)
            km = np.array(kc)
        return np.sort(km)

    for k in range(3,8):
        print(k)
        kmers = generatek(k)
        print(kmers)
        kfeat = np.zeros((len(seqs), len(kmers)), dtype = np.float16)
        for s, seq in enumerate(ids):
            seq = seq.replace('seq_', '')
            kset = []
            for l in range(len(seq)-k+1):
                kset.append(seq[l:l+k])
            uk, un = np.unique(kset, return_counts = True)
            kfeat[s, np.isin(kmers, uk)] = un
        kmers = np.array(['bar_'+k for k in kmers])
        np.savez_compressed('Barcodekmer'+str(k)+'_input.npz', names = ids, features = kfeat, featuretypes = kmers)
        print('Barcodekmer'+str(k)+'_input.npz')

 

