
# coding: utf-8


import os,sys
from scipy.stats.stats import pearsonr
from scipy.cluster.vq import kmeans, vq
import numpy as np

pairwise = []

def average_match(e, l, ignore_self=False):
    """The average match between element e and list l"""
    # filter out ourself as this will give us a nice 1.0 value
    # the problem with this is if we are a list of 1, anything will have a better value!
    if ignore_self:
        l = [x for x in l if x != e]
    if len(l) == 0:
        return 0
    values = [pairwise[e][x] for x in l]
    return 1.0*sum(values)/len(values)



#data_file = "output.contigs2reads.contigs.short.txt"
try :
    data_file = sys.argv[1]
except:
    sys.exit(sys.argv[0] + " <contigs file>")

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

# measure all correlations
for i in range(len(ids)):
    for j in range(i+1, len(ids)):
        pearson, p = pearsonr(data[ids[i]], data[ids[j]])
        if pearson == np.nan:
            continue
        pairwise[i][j] = pearson



# generate an initial estimate of the clusters by using a clustering 
# algorithm to group things with pearson > 0.9

cluster = {}
members = []
cluster_count=0

# go through and find all matches > 90. Add them to sets
# we only add a member if adding it would not reduce the average percent
# id of the cluster below 80. This is more relaxed than the overall
# percent id
for i in range(len(ids)):
    # we need a cluster for this item if we don't have one
    if i not in cluster:
        cluster[i] = cluster_count
        members.append([i])
        cluster_count += 1
    for j in range(i+1, len(ids)):
        if j in cluster:
            continue
        if pairwise[i][j] > 0.9:
            cluster[j] = cluster[i]
            members[cluster[i]].append(j)
            # has this reduced our average too far?
            if average_match(j, members[cluster[j]]) > 0.8:
                cluster.pop(j)
                members[cluster[i]].pop()


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
