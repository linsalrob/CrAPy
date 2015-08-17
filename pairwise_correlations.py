

# coding: utf-8


"""
generate all pairwise correlations and print a list of contig1, contig2, correlation score
"""

import os,sys
from scipy.stats.stats import pearsonr
import numpy as np
from multiprocessing import Pool, Lock


NUM_THREADS=6



def pairwise_correlation(data, ids, i, j):
        try:
            pearson, p = pearsonr(data[ids[i]], data[ids[j]])
            if pearson == np.nan:
                pearson = 0
                p = 0
        except Exception as e:
            sys.stderr.write("There was an error when i: " + str(i) + " and j " + str(j) + " messsage: " + str(e) + "\n")
        return [ids[i], ids[j], pearson, p]





if __name__ == '__main__':
    #data_file = "test_data/output.contigs2reads.contigs.100.txt"
    try:
        data_file = sys.argv[1]
        noc = float(sys.argv[2])
    except:
        sys.exit(sys.argv[0] + " <data file e.g. test_data/output.contigs2reads.contigs.100.txt> <minimum number of occurences (try 3)>")

    if noc < 3:
        sys.stderr.write("Warning: you choose reads that have at least " + str(noc) + " occurrences, but this might be too low. Continuing anyway\n")


    data = {}
    with open(data_file, 'r') as fin:
        for l in fin:
            p=l.strip().split("\t")
            did = p.pop(0)
            p = map(float, p)
            # only keep those that have > 3 zeros
            if len(p) - p.count(0) <= noc:
                continue
            data[did] = p

    # now we need to find all pairwise correlation scores
    ids = data.keys()

    # initial cluster
    pool = Pool(NUM_THREADS)

    # measure all correlations
    results = [pool.apply_async(pairwise_correlation, [data, ids, i, j]) for i in range(len(ids)) for j in range(i+1, len(ids))]
    results = [p.get() for p in results]
    for r in results:
        print("\t".join(map(str, r)))


