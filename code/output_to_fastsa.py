import numpy as np
import sys, os

obj = open(sys.argv[1], 'r').readlines()
out = open(os.path.splitext(sys.argv[1])[0]+'.fasta', 'w')
for l, line in enumerate(obj):
    if line[0] == '>':
        out.write(line+obj[l+1]+obj[l+2].split()[0].replace('+','.')+'\n')

