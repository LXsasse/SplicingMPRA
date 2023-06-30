import numpy as np
import sys, os


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
        out = line.strip().split('\t').index(outputfeat)
        sc = line.strip().split('\t').index('sequence')
    else:
        line = line.strip().split('\t')
        if len(line[0]) > 10:
            ids.append('seq_'+line[0])
            ef.append(line[out])
            seqs.append(line[sc])

ids = np.array(ids)
ef = np.array(ef)
seqs = np.array(seqs)

sort = np.argsort(ids)
ids, ef, seqs = ids[sort], ef[sort], seqs[sort]
nts = np.array(list('ACGT'))
onehot = []
for seq in seqs:
    onehot.append(np.array(list(seq))[:,None] == nts)
onehot = np.array(onehot)

np.savetxt('N50_meanefficacy.csv', np.array([ids, ef]).astype(str).T, fmt = '%s', delimiter = ',')
np.savez_compressed('N50_seqs.npz', genenames = ids, seqfeatures = onehot, featurenames = nts)
 

