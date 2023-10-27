*Pseudocode for BFR*

take k centroids (c1, c2, ..., ck) in the dataset space

init k discarder sets (D1, D2, ..., Dk) {number of points, vector of the sum of coordinates, vector of the sum of squares of coordinates} 
to {0, [0* d], [0* d]}

start streaming the dataset

while (stream is not finished){
    take a new point p from the stream

    if (p is closer then threshold to some ci){
        update the discarder set Di with p
    }else if (p is closer then threshold to some compressed set Cj){
        update the compressed set Cj with p
    }else{
        add p to the retained set R (R = R U {p})
    }
}

update centroids (c1, c2, ..., ck) with the statistics of each discarder set (D1, D2, ..., Dk)

start merging compressed set till convergence:
    use statistics to decide if two sets should be merged (e.g. if with new centroid variance is smaller than threshold)

start merge retained set with cluster using malahanobis distance:
    use statistics to decide if a point should be merged with a cluster (e.g. if malahanobis distance between point and cluster centroid is smaller than threshold)

return the final clusters centroids