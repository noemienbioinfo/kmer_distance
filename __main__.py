from cProfile import label


def main():
    import distance as dist
    from skbio import DistanceMatrix
    from skbio.tree import nj

    #whole k-mers dictionaries
    R6_a = dist.kmers_from_seq("R6.fa", 14)
    TIGR_a = dist.kmers_from_seq("TIGR4.fa", 14)
    draft82_a = dist.kmers_from_seq("14412_3#82.contigs_velvet.fa", 14)
    draft84_a = dist.kmers_from_seq("14412_3#84.contigs_velvet.fa", 14)

    #simple Jaccard distances
    jacc1 = dist.jaccard_dist(R6_a, TIGR_a)
    jacc2 = dist.jaccard_dist(R6_a, draft82_a)
    jacc3 = dist.jaccard_dist(R6_a, draft84_a)
    jacc4 = dist.jaccard_dist(TIGR_a, draft82_a)
    jacc5 = dist.jaccard_dist(TIGR_a, draft84_a)
    jacc6 = dist.jaccard_dist(draft82_a, draft84_a)

    print("Simple Jaccard distances with whole dictionaries :\n")
    print(f"R6 vs TIGR4 : {jacc1}")
    print(f"R6 vs draft82 : {jacc2}")
    print(f"R6 vs draft84 : {jacc3}")
    print(f"TIGR4 vs draft82 : {jacc4}")
    print(f"TIGR4 vs draft84 : {jacc5}")
    print(f"draft82 vs draft84 : {jacc6}")

    #loading of text files
    R6_hash = dist.file2list("R6.txt")
    TIGR_hash = dist.file2list("TIGR4.txt")
    draft82_hash = dist.file2list("draft82.txt")
    draft84_hash = dist.file2list("draft84.txt")

    R6vsTIGR4 = dist.minhash_jaccard(R6_hash, TIGR_hash)
    R6vsdraft82 = dist.minhash_jaccard(R6_hash, draft82_hash)
    R6vsdraft84 = dist.minhash_jaccard(R6_hash, draft84_hash)
    TIGR4vsdraft82 = dist.minhash_jaccard(TIGR_hash, draft82_hash)
    TIGR4vsdraft84 = dist.minhash_jaccard(TIGR_hash, draft84_hash)
    draft82vsdraft84 = dist.minhash_jaccard(draft82_hash, draft84_hash)

    dist_list = [R6vsTIGR4, R6vsdraft82, R6vsdraft84, TIGR4vsdraft82, TIGR4vsdraft84, draft82vsdraft84]

    #minhash Jaccard distances. Use of text files
    print("\nMinHash Jaccard distances :\n")
    print(f"R6 vs TIGR4 : {R6vsTIGR4}")
    print(f"R6 vs draft82 : {R6vsdraft82}")
    print(f"R6 vs draft84 : {R6vsdraft84}")
    print(f"TIGR4 vs draft82 : {TIGR4vsdraft82}")
    print(f"TIGR4 vs draft84 : {TIGR4vsdraft84}")
    print(f"draft82 vs draft84 : {draft82vsdraft84}")

    #neighbor joining tree. Calculation of the distance matrix then draw tree
    ids = ['R6', 'TIGR4', 'draft82', 'draft84']
    dm = DistanceMatrix(dist_list, ids)

    tree = nj(dm)
    print()
    print(tree.ascii_art())



def multiple(kmers=14, sketch_min=100, sketch_max=1000, step=100):
    import distance as dist
    import matplotlib.pyplot as plt

    #calculate distances with Minhash calc from raw data and append distances into
    #lists in chosen range
    x = [0, ]
    y_R6_TIGR = [0, ]
    y_R6_D82 = [0, ]
    y_R6_D84 = [0, ]
    y_TIGR_D82 = [0, ]
    y_TIGR_D84 = [0, ]
    y_D82_D84 = [0, ]
    R6_a = dist.kmers_from_seq("R6.fa", kmers)
    TIGR_a = dist.kmers_from_seq("TIGR4.fa", kmers)
    draft82_a = dist.kmers_from_seq("14412_3#82.contigs_velvet.fa", kmers)
    draft84_a = dist.kmers_from_seq("14412_3#84.contigs_velvet.fa", kmers)
    for i in range(sketch_min, sketch_max, step):
        x.append(i)
        hashed_R6 = dist.convert_hash(R6_a, i)
        hashed_TIGR = dist.convert_hash(TIGR_a, i)
        hashed_D82 = dist.convert_hash(draft82_a, i)
        hashed_D84 = dist.convert_hash(draft84_a, i)
        Hdist_R6_TIGR = dist.minhash_jaccard(hashed_R6, hashed_TIGR)
        Hdist_R6_D82 = dist.minhash_jaccard(hashed_R6, hashed_D82)
        Hdist_R6_D84 = dist.minhash_jaccard(hashed_R6, hashed_D84)
        Hdist_TIGR_D82 = dist.minhash_jaccard(hashed_TIGR, hashed_D82)
        Hdist_TIGR_D84 = dist.minhash_jaccard(hashed_TIGR, hashed_D84)
        Hdist_D82_D84 = dist.minhash_jaccard(hashed_D82, hashed_D84)
        y_R6_TIGR.append(Hdist_R6_TIGR)
        y_R6_D82.append(Hdist_R6_D82)
        y_R6_D84.append(Hdist_R6_D84)
        y_TIGR_D82.append(Hdist_TIGR_D82)
        y_TIGR_D84.append(Hdist_TIGR_D84)
        y_D82_D84.append(Hdist_D82_D84)
       
    #from x list and y lists, plot graph then draw it
    fig, ax = plt.subplots()

    ax.plot(x, y_R6_TIGR, label='R6vsTIGR')
    ax.plot(x, y_R6_D82, label='R6vsD82')
    ax.plot(x, y_R6_D84, label='R6vsD84')
    ax.plot(x, y_TIGR_D82, label='TIGRvsD82')
    ax.plot(x, y_TIGR_D84, label='TIGRvsD84')
    ax.plot(x, y_D82_D84, label='D82vsD84')
    plt.legend()
    plt.show()
        


if __name__ == "__main__":
    main()
    # multiple()