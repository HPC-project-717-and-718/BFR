*Pseudocode for BFR*

take k centroids (c1, c2, ..., ck) in the dataset space

init k discarder sets (DS1, DS2, ..., DSk) {number of points, vector of the sum of coordinates, vector of the sum of squares of coordinates} 
to {0, [0* d], [0* d]}

start streaming the dataset D (size |D|), with each iteration take a block B of size |B| < |D|. B should fill the primary memory (RAM).

while (dataset not read completely){
    for all points p in B{
        if (p is closer than threshold to some ci or second primary compression criteria is found){ 
            update the discarder set Di with p
        }else {
            add p to the retained set R (R = R U {p})
        }
    }

    update centroids (c1, c2, ..., ck) with the statistics of each discarder set (D1, D2, ..., Dk)

    cluster retained set R with classical K-Means, creating k2 clusters (k2>K):
        if (new cluster k2[i] fits criteria (has more than one point OR fits a tightness criteria)){
            add k2[i] to list of compressed sets
        }else {
            add k2[i] back to retained set R
        }

    start merging compressed sets till convergence using hierarchical agglomerative clustering:
        use statistics to decide if two sets should be merged (e.g. if with new centroid variance is smaller than threshold)

    take a new block B
}

decide what to do with remaining compressed sets and retained set:
    treat them as outliers.

return the final clusters centroids