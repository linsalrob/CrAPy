"""
Convert pairwise correlations to clusters
"""

import sys
import os



def average_match(pairwise, e, l, ignore_self=False):
    """The average match between element e and list l"""
    # filter out ourself as this will give us a nice 1.0 value
    # the problem with this is if we are a list of 1, anything will have a better value!
    if ignore_self:
        l = [x for x in l if x != e]
    if len(l) == 0:
        return 0
    values = [pairwise[e][x] for x in l]
    return 1.0*sum(values)/len(values)






if __name__ == "__main__":
    try:
        pairwiseF = sys.argv[1]
    except:
        sys.exit(sys.argv[0] + " <pairwise correlation score (output from pairwise.py)>")


    pairwise = {}
    with open(pairwiseF, 'r') as pin:
        for l in pin:
            p=l.strip().split('\t')
            if p[0] not in pairwise:
                pairwise[p[0]]={}
            if p[1] not in pairwise:
                pairwise[p[1]]={}
            pairwise[p[0]][p[1]]=float(p[2])
            pairwise[p[1]][p[0]]=float(p[2])
    allids = list(pairwise.keys())
    allids.sort()


    # generate the clusters, with an initial greedy clustering 
    # at 0.90
    clustercount = 0
    cluster = {}
    members = {}
    for i in range(len(allids)):
        idi = allids[i]
        if i not in cluster:
            cluster[i] = clustercount
            members[clustercount]={i}
            clustercount+=1
        for j in range(i+1, len(allids)):
            if j in cluster:
                continue
            idj = allids[j]
            if pairwise[idi][idj] > 0.9:
                cluster[j]=cluster[i]
                members[cluster[i]].add(j)

    # for any element in the list, I want to find if it matches another cluster better than its own cluster


    for i in range(len(allids)):
        best_score = 0
        best_id = None
        iclust = cluster[i]
        memberset = members[iclust]

        # compare to the other clusters
        for c in len(members):
            if (c == icluster):
                continue
            score = average_match(pairwise, i, members[c], True)

            if score > best_score and score > 0:
                best_score = score
                best_id = c

        memberofscore = average_match(pairwise, i, members[iclust], True)

        if best_score > memberscore:
            cluster[i] = best_id
            members[iclust].pop(i)
            members[best_id].add(i)


    for i in cluster:
        print("\t".join(map(str, [i, cluster[i]])))
                        


