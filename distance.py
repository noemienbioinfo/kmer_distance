from Bio import SeqIO
from mmh3 import hash

def kmers_into_dict(file):
    list_seq = []
    for record in SeqIO.parse(file, "fasta"):
        list_seq.append(record.seq)
    mers = dict()
    for seq in list_seq:
        length = len(seq)-13
        for i in range(length):
            k = hash(str(seq[i:i+14]))
            if k in mers:
                mers[k] += 1
            else:
                mers[k] = 1
    return mers

# print(kmers_into_dict("R6.fa"))
# print(kmers_into_dict("TIGR4.fa"))