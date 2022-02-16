from cgitb import text
from Bio import SeqIO
from mmh3 import hash

#STEP1
def kmers_into_dict(file):
    #make a list of sequences in the file (i.e. draft sequences with multiple seq)
    list_seq = []
    for record in SeqIO.parse(file, "fasta"):
        list_seq.append(record.seq)

    #make a dictionary with the hashed 14-mers associated with the number it's counted (the number of counts will not be used)
    mers = dict()
    for seq in list_seq:
        length = len(seq)-13
        for i in range(length):
            k = hash(str(seq[i:i+14])) #format(str, 'X').lower()
            if k in mers:
                mers[k] += 1
            else:
                mers[k] = 1
    #mers = dict(hashed(seq) : nbr_of_occurences)
    return mers

# print(kmers_into_dict("R6.fa"))
# print(kmers_into_dict("TIGR4.fa"))

#STEP2
def jaccard(seq1, seq2):
    #make a dictionary of commun mers between seq1 and seq2
    #it's not mandatory to add seq1[mers] and seq2[mers]
    commun_mers = dict()
    for mers in seq1:
        if mers in seq2:
            commun_mers[mers] = seq1[mers] + seq2[mers]

    #count of union unique mers ("if" is necessary to avoid duplicated mers)
    union_mers = 0
    for mers in seq1:
        union_mers += 1
    for mers in seq2:
        if mers not in seq1:
            union_mers += 1

    #count of inter mers which is the number of mers in the commun_mers dict
    inter_mers = 0
    for mers in commun_mers:
        inter_mers += 1
    
    #1 - J(A, B)
    return (str(inter_mers), str(union_mers))

# R6 = kmers_into_dict("R6.fa")
# TIGR = kmers_into_dict("TIGR4.fa")
# print(R6)
# print(jaccard(R6, TIGR))

# hash_run1 = kmers_into_dict("R6.fa")
# hash_run2 = kmers_into_dict("R6.fa")

# if hash_run1 == hash_run2:
#     print("ok")
# else:
#     print("NOPE")

#STEP3
def sketch(seq):
    #put unique mers in a liste before sorting and get the first thousand of mers
    list_mers = []
    for elem in seq:
        list_mers.append(elem)
    list_mers.sort()
    list_mers = list_mers[0:1000]

    #turn liste of first thousand mers into dict() for Jaccard fonction (which uses dict())
    #the value of keys has no importance because we don't use the number of occurences of each mers
    dico_mers = dict()
    for elem in list_mers:
        dico_mers[elem] = 1

    #dico of first thousand of mers (sorted by hashes), ready to use for jaccard(seq1, seq2)
    return dico_mers

# sketch_r6 = sketch(R6)
# sketch_tigr = sketch(TIGR)
# J_all = jaccard(R6, TIGR)
# J_sketch = jaccard(sketch_r6, sketch_tigr)
# print("Jaccard R6/TIGR4 =", J_all[0] + "/" + J_all[1], "\n"
#     "Jaccard of sketches R6/TIGR4 =", J_sketch[0] + "/" + J_sketch[1]
#     )
# print(sketch_r6)

list_outputs = ["R6.txt", "TIGR4.txt", "14412_3#82.txt", "14412_3#84.txt"]
list_fasta = ["R6.fa", "TIGR4.fa", "14412_3#82.contigs_velvet.fa", "14412_3#84.contigs_velvet.fa"]
def edit(list_outputs, list_inputs):
    for files in range(len(list_outputs)):
        with open(list_outputs[files], "w") as opf:
            text_to_write = sketch(kmers_into_dict(list_inputs[files]))
            print(str(text_to_write), file=opf)

# edit(list_outputs, list_fasta)

with open("R6.txt", "r") as file1:
    with open("TIGR4.txt", "r") as file2:
        file1 = file1.read()
        file2 = file2.read()
        for k, v in file1:
            v = int(v)
        for k, v in file2:
            v = int(v)
        Jaccard = jaccard(file1, file2)
        print("Jaccard of sketches R6/TIGR4 from sketches files =", Jaccard[0]\
            + "/" + Jaccard[1])
        
#Test for draft sequences. Question : are the sequences concatenated ?
# test1 = print(kmers_into_dict("test.fa"))
#answer : if multiple sequences in FASTA file, each seq is independant