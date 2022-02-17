from Bio import SeqIO
import mmh3

#STEP1
def kmers_from_seq(file, k):
    #make a list of sequences in the file (i.e. draft sequences with multiple seq)
    list_seq = []
    for record in SeqIO.parse(file, "fasta"):
        list_seq.append(record.seq)

    #make a dictionary with the hashed 14-mers associated with the number it's counted (the number of counts will not be used)
    mers = dict()
    for seq in list_seq:
        length = len(seq)-k-1
        for i in range(length):
            kmer = seq[i:i+k]
            if kmer in mers:
                mers[kmer] += 1
            else:
                mers[kmer] = 1

    #mers = dict(14-mers : nbr_of_occurences)
    return mers

#STEP2
def jaccard_dist(seq1, seq2):
    #make a dictionary of commun mers between seq1 and seq2
    #it's not mandatory to add seq1[mers] and seq2[mers]
    commun_mers = dict()
    for mers in seq1:
        if mers in seq2:
            commun_mers[mers] = seq1[mers] + seq2[mers]

    #count of union unique mers ("if" is necessary to avoid duplicated mers)
    union = len(seq1)
    for mers in seq2:
        if mers not in seq1:
            union += 1

    #count of inter mers which is the number of mers in the commun_mers dict
    inter = len(commun_mers)
    
    #1 - J(A, B)
    return 1-inter/union

#STEP3
#generate a list of all unique k-mers, then sort it and return the firsts thousands
def convert_hash(dict_kmers, size):
    list_kmers = []
    for mers in dict_kmers:
        hashed = mmh3.hash(str(mers), signed = False)
        list_kmers.append(hashed)
    list_kmers.sort()

    return list_kmers[0:size]

#this function is like jaccard_dist function but uses lists (not dict)
def minhash_jaccard(list_kmers1, list_kmers2):
    commun_mers = []
    for mers in list_kmers1:
        if mers in list_kmers2:
            commun_mers.append(mers)
    
    union = len(list_kmers1)
    for mers in list_kmers2:
        if mers not in list_kmers1:
            union += 1
    
    inter = len(commun_mers)

    return 1-inter/union

#STEP4
#function used to edit lists of hashed kmers
def edit(list_outputs, list_inputs):
    for files in range(len(list_outputs)):
        with open(list_outputs[files], "w") as opf:
            sketch = convert_hash(kmers_from_seq(list_inputs[files]))
            for hashed_kmer in sketch:
                opf.write(str(hashed_kmer) + "\n")

#function to read text files and translate into lists
def file2list(file):
    kmers_list = []
    with open (file, "r") as ipf:
        for line in ipf:
            kmers_list.append(line.rstrip())

    return kmers_list

#make distance matrix used for drawing tree
def make_matrix(list_of_scores):
    matrix = []
    for i in range(len(list_of_scores)):
        line_matrix = []
        for j in range(len(list_of_scores)):
            if list_of_scores[i] >= list_of_scores[j]:
                line_matrix.append(list_of_scores[i]-list_of_scores[j])
            else:
                line_matrix.append(list_of_scores[j]-list_of_scores[i])
        matrix.append(line_matrix)
    return matrix