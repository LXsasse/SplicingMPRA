import numpy as np
import sys, os
import glob
def readps(psfile):
    obj = open(psfile, 'r').readlines()
    seq = None
    readseq = False
    for l, line in enumerate(obj):
        if readseq:
            if line[:3] == ') }':
                readseq = False
                matrix = np.zeros((len(seq), len(seq)))
            else:
                seq += line.strip().strip("\\")

        if line[:9] == '/sequence':
            readseq = True
            seq = ''

        if seq is not None and readseq == False and line.strip()[-4:] == 'ubox':
            line = line.strip().split()
            matrix[int(line[0])-1, int(line[1])-1] = float(line[2])**2
            matrix[int(line[1])-1, int(line[0])-1] = float(line[2])**2
    return seq, matrix



folder = sys.argv[1]
psfiles = np.sort(glob.glob(folder+'*dp.ps'))

makemat = False
if '--generate_matrix_rep' in sys.argv:
    matrices = []
    makemat = True
seqs = []
acc = []
names = []
nucs = np.array(list('ACGT'))
for f, fi in enumerate(psfiles):
    if f % 100 == 0:
        print(f)
    seq, matrix = readps(fi)
    name = fi.split('_dp.ps')[0]
    if '/' in name:
        name = name.split('/')[1]
    names.append(name)
    seq = seq.replace('U', 'T')
    seq = np.array(list(seq))[:,None] == nucs
    seqs.append(seq)
    acc.append(np.around(np.sum(matrix, axis = -1),3).astype(np.float16))
    if makemat:
        matrices.append(np.around(matrix,3).astype(np.float16))

seqs = np.array(seqs).astype(np.float32)
if makemat:
    matrices = np.array(matrices).astype(np.float16)
    matrix_features = np.append(seqs, matrices, axis = -1)
    np.savez_compressed('Matrix_representation.npz', seqfeatures = matrix_features, genenames = np.array(names), featurenames = np.append(nucs, np.arange(len(seqs[0])).astype(str)))

access_features = np.append(seqs, access[...,None], axis = -1)

np.savez_compressed('Access_representation.npz', seqfeatures = access_features, genenames = np.array(names), featurenames = np.append(nucs, ['access']))
