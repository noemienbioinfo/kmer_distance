import distance

#######################################################################
##########   prints to test the production of dictionaries   ##########
##########              by kmers_from_seq function           ##########
#######################################################################

# print(distance.kmers_from_seq("R6.fa"))
# print(distance.kmers_from_seq("TIGR4.fa"))

#######################################################################
##########      prints to test the distance calculation      ##########
#######################################################################

# R6 = distance.kmers_from_seq("R6.fa")
# TIGR = distance.kmers_from_seq("TIGR4.fa")
# draft82 = distance.kmers_from_seq("14412_3#82.contigs_velvet.fa")
# draft84 = distance.kmers_from_seq("14412_3#84.contigs_velvet.fa")
# print(distance.jaccard_dist(R6, TIGR))
# print(distance.jaccard_dist(draft82, draft84))

#######################################################################
##########     verification of same hash for same sample     ##########
#######################################################################

# hash_run1 = distance.kmers_into_dict("R6.fa")
# hash_run2 = distance.kmers_into_dict("R6.fa")

# if hash_run1 == hash_run2:
#     print("ok")
# else:
#     print("NOPE")

#######################################################################
##########    prints to test the minhash_jaccard function    ##########
##########         and compare to jaccard_dist function      ##########
#######################################################################

# sketch_r6 = distance.minhash_jaccard(R6)
# sketch_tigr = distance.minhash_jaccard(TIGR)
# J_all = distance.jaccard_dist(R6, TIGR)
# J_sketch = distance.jaccard_dist(sketch_r6, sketch_tigr)
# print("Jaccard R6/TIGR4 =", J_all[0] + "/" + J_all[1], "\n"
#     "Jaccard of sketches R6/TIGR4 =", J_sketch[0] + "/" + J_sketch[1]
#     )
# print(sketch_r6)

#######################################################################
##########         edition of txt files from sketches        ##########
#######################################################################

# list_outputs = ["R6.txt", "TIGR4.txt", "draft82.txt", "draft84.txt"]
# list_fasta = ["R6.fa", "TIGR4.fa", "14412_3#82.contigs_velvet.fa", "14412_3#84.contigs_velvet.fa"]
# distance.edit(list_outputs, list_fasta)

#######################################################################
##########        test minhash_jaccard from .txt files       ##########
#######################################################################

file1 = distance.file2list("R6.txt")
file2 = distance.file2list("TIGR4.txt")
print(distance.minhash_jaccard(file1, file2))

#######################################################################
##########             test for the drafts files             ##########
#######################################################################

#Question : are the sequences concatenated ?
# test1 = print(kmers_into_dict("test.fa"))
#answer : if multiple sequences in FASTA file, each seq is independant