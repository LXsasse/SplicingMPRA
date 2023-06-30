import numpy as np
import sys, os

nucs = np.array(list('ACGT'))

sfasta = open(sys.argv[1], 'r').readlines()
names = []
seqs = []
brackets = []
for l, line in enumerate(sfasta):
    if line[0] == '>':
        names.append(line[1:].strip())
        seq = sfasta[l+1].strip()
        seq = seq.replace('U', 'T')
        seq = np.array(list(seq))[:,None] == nucs
        seqs.append(seq)
        brackets.append(sfasta[l+2].strip())

mfasta = open(sys.argv[2], 'r').readlines()
mnames = []
struc = []
se = np.array(list('fshmi'))
i = 0
for l, line in enumerate(mfasta):
    if line[:1] == '>':
        mnames.append(line[1:].strip())
        if mnames[-1] == names[i]:
            st = mfasta[l+1].strip().replace('t', 'f')
            st = np.array(list(st))[:,None] == se
            struc.append(st)
            i += 1
        else:
            print('BRACKETS dont match', i, l)
            print(mnames[-1])
            print(names[i])
            sys.exit()

seqs = np.array(seqs).astype(np.float32)
struc = np.array(struc).astype(np.float32)
print(np.shape(seqs), np.shape(struc))
struc_features = np.append(seqs, struc, axis = -1)
print(np.shape(struc_features))

np.savez_compressed('SeqStructure_representation.npz', seqfeatures = struc_features, genenames = np.array(names), featurenames = np.append(nucs, se))
