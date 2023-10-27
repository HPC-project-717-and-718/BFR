import numpy as np
import sys
import os

output_file = sys.argv[1]
dimensions = sys.argv[2]

if not int(dimensions) >= 2:
    sys.exit("Dimensions must be 2 or higher.")

number_of_points = int(dimensions) * 20
number_of_gaussians = 10
minimum_mean, maximum_mean = -5, 5
minimum_covariance, maximum_covariance = 0.7, 1.5

output = open(output_file, "w")
# write number of points and number of dimensions
output.write(str(number_of_points) + "\t" + str(dimensions) + "\n")

# Step 1. Generate number_of_gaussian means and covariance matrixes
# generate number_of_gaussians means
means = np.random.uniform(minimum_mean, maximum_mean, (number_of_gaussians, int(dimensions)))

# generate number_of_gaussians covariance matrixes with size dimensions x dimensions and matrixes are diagonal
covariances_diagonal = np.zeros((number_of_gaussians, int(dimensions), int(dimensions)))
covariances = np.random.uniform(minimum_covariance, maximum_covariance, (number_of_gaussians, int(dimensions)))
for i in range(number_of_gaussians):
    covariances_diagonal[i] = np.diag(covariances[i])

# Step 2. Generate number_of_points points
points = np.zeros((number_of_points, int(dimensions)))

# give weight to each gaussian
weights = np.random.uniform(0, 1, number_of_gaussians)
weights = weights / np.sum(weights)

# initialize weights count
weights_count = np.zeros(number_of_gaussians)

for i in range(number_of_points):
    # choose a gaussian by weight
    gaussian = np.random.choice(number_of_gaussians, p=weights)

    #count how many times each gaussian is chosen
    weights_count[gaussian] += 1

    # generate a point from the gaussian
    points[i] = np.random.multivariate_normal(means[gaussian], covariances_diagonal[gaussian])

    # print point in output file
    for j in range(int(dimensions)):
        output.write(str(points[i][j]) + "\t")
    output.write("\n")

# print weights count
# print(weights_count)