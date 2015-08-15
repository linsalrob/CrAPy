
# coding: utf-8


import os,sys
from scipy.stats.stats import pearsonr
import numpy as np
from multiprocessing import Pool, Lock


NUM_THREADS=6

pairwise = []
lock = None
cluster_count = 0
cluster = {}

def init(l, c, c_c):
    global lock 
    lock = l
    global cluster
    cluster = c
    global cluster_count 
    cluster_count = c_c


def average_match(e, l, ignore_self=False):
    """The average match between element e and list l"""
    # filter out ourself as this will give us a nice 1.0 value
    # the problem with this is if we are a list of 1, anything will have a better value!
    global pairwise
    if ignore_self:
        l = [x for x in l if x != e]
    if len(l) == 0:
        return 0
    values = [pairwise[e][x] for x in l]
    return 1.0*sum(values)/len(values)



def pairwise_correlation(data, ids, i, j):
        try:
            pearson, p = pearsonr(data[ids[i]], data[ids[j]])
            if pearson == np.nan:
                pearson = 0
        except Exception as e:
            sys.stderr.write("There was an error when i: " + str(i) + " and j " + str(j) + " messsage: " + str(e) + "\n")
        if j == 500:
            print(str(i) + " : " + str(j))
        return [i, j, pearson]


def generate_clusters(pairwise, i, j):

    # go through and find all matches > 90. Add them to sets
    # we only add a member if adding it would not reduce the average percent
    # id of the cluster below 80. This is more relaxed than the overall
    # percent id

    global cluster_count
    global lock
    global cluster

    if j in cluster:
        return

    if i not in cluster:
        lock.acquire()
        try:
            cluster[i] = cluster_count
            cluster_count += 1
        finally:
            lock.release()

    if pairwise[i][j] > 0.9:
        lock.acquire()
        try:
            cluster[j] = cluster[i]
        finally:
            lock.release()






if __name__ == '__main__':
    data_file = "test_data/output.contigs2reads.contigs.100.txt"

    data = {}
    with open(data_file, 'r') as fin:
        for l in fin:
            p=l.strip().split("\t")
            did = p.pop(0)
            p = map(float, p)
            data[did] = p

    # now we need to find all pairwise correlation scores
    ids = data.keys()

    for i in ids:
        temp = []
        for j in ids:
            temp.append(0)
        pairwise.append(temp)

    # set ourselves to 1
    for i in range(len(ids)):
        pairwise[i][i] = 1

    # how many processes to use
    l = Lock()
    pool = Pool(NUM_THREADS, initializer=init, initargs=(l, {}, 0, ))
    # measure all correlations
    results = [pool.apply_async(pairwise_correlation, [data, ids, i, j]) for i in range(len(ids)) for j in range(i+1, len(ids))]
    results = [p.get() for p in results]

    for r in results:
        pairwise[r[0]][r[1]]=r[2]
        pairwise[r[1]][r[0]]=r[2]


    # generate an initial estimate of the clusters by using a clustering 
    # algorithm to group things with pearson > 0.9

    
    results = [pool.apply_async(generate_clusters, [pairwise, i, j]) for i in range(len(ids)) for j in range(i+1, len(ids))]
    results = [p.get() for p in results]

    sys.stderr.write("After assigning the clusters there are " + str(cluster_count) + " clusters or len: " + str(len(cluster)) + "\n")

    members = []
    for i in range(cluster_count):
        members.append([])
    for i in cluster:
        members[cluster[i]].append(i)

    # for any element in the list, I want to find if it matches another cluster better than its own cluster


    for i in range(len(ids)):
        best_score = 0
        best_id = None
        for j in range(len(members)):
            score = average_match(i, members[j], True)
            if score > best_score and score > 0:
                best_score = score
                best_id = j
        if cluster[i] != best_id:
            curr_score = average_match(i, members[cluster[i]], True)
            #print(str(i) + " current: " + str(cluster[i]) + " best: " + str(best_id) + " current score: " + str(curr_score) + " best score: " + str(best_score))
            cluster[i] = best_id

    for i in range(len(ids)):
        if i not in cluster or not cluster[i]:
            cluster[i] = cluster_count
            cluster_count += 1
        print("\t".join(map(str, [ids[i], cluster[i]])))
