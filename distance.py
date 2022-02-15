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
            k = format(hash(str(seq[i:i+14])), 'X').lower()
            if k in mers:
                mers[k] += 1
            else:
                mers[k] = 1
    return mers

# print(kmers_into_dict("R6.fa"))
# print(kmers_into_dict("TIGR4.fa"))

def jaccard(seq1, seq2):
    commun_mers = dict()
    for mers in seq1:
        if mers in seq2:
            commun_mers[mers] = seq1[mers] + seq2[mers]

    union_mers = 0
    for mers in seq1:
        union_mers += 1
    for mers in seq2:
        if mers not in seq1:
            union_mers += 1

    inter_mers = 0
    for mers in commun_mers:
        inter_mers += 1
    
    return 1-(inter_mers/union_mers)

R6 = kmers_into_dict("R6.fa")
# print(R6)
TIGR = kmers_into_dict("TIGR4.fa")
print(jaccard(R6, TIGR))

# hash_run1 = kmers_into_dict("R6.fa")
# hash_run2 = kmers_into_dict("R6.fa")

# if hash_run1 == hash_run2:
#     print("ok")
# else:
#     print("NOPE")