import numpy as np
import sys, os
import matplotlib.pyplot as plt

def checkint(i):
    try: 
        i = int(i)
        return True
    except:
        return False

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


seq, matrix = readps(sys.argv[1])

fig = plt.figure(figsize = (len(seq)*0.2, len(seq)*0.2 *1.06), dpi = 20)
ax = fig.add_subplot(211)
ax.set_position([0.1,0.1,0.7,0.7])
ax.imshow(matrix, cmap = 'Greys', vmin = 0, vmax = 1, aspect = 'auto')
ax.set_xticks(np.arange(len(seq)))
ax.set_xticklabels(list(seq))
ax.set_yticks(np.arange(len(seq)))
ax.set_yticklabels(list(seq))

ax2 = fig.add_subplot(212)
ax2.set_position([0.1,0.81,0.7,0.05])
ax2.bar(np.arange(len(seq)), np.sum(matrix,axis = 1), color='k')
ax2.set_xlim([0,len(seq)])

fig.savefig(os.path.splitext(sys.argv[1])[0]+'.jpg', dpi = 100, bbox_inches = 'tight')

plt.show()
