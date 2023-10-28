import matplotlib.pyplot as plt
import numpy as np
import sys

output_file = sys.argv[1]
dimensions_str = sys.argv[2]
plotname = sys.argv[3]
dimensions = int(dimensions_str)

output = open(output_file, "r")
firstline = output.readline()
points = output.readlines()

# dimensions_axis holds the points' coordinates for each axis
dimension_axis = []
for i in range(dimensions):
    dimension_axis.append([])


# transform points to float
for i in range(len(points)):
    points[i] = points[i].split("\t")
    # remove element containing \n
    points[i].pop()
    for j in range(dimensions):
        points[i][j] = float(points[i][j])
        dimension_axis[j].append(points[i][j])
        
# plot the points
# plt.figure()
# plt.scatter(dimension_axis[0], dimension_axis[1], marker='o', color='b', s=0.5, label='Data Points')
# plt.title('Data Plot')
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.show()
# plt.savefig(plotname)

# plot the points in 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(dimension_axis[0], dimension_axis[1], dimension_axis[2], c='r', s=3, marker='o')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()
plt.savefig(plotname)