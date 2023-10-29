import matplotlib.pyplot as plt
import sys
import matplotlib.patches as patches
import math

output_file = sys.argv[1]
dimensions_str = sys.argv[2]
plotname = sys.argv[3]
dimensions = int(dimensions_str)

output = open(output_file, "r")
firstline = output.readline()
points = output.readlines()

# dimensions_axis holds the points' coordinates for each axis
dimension_axis = [[] for _ in range(dimensions)]

# transform points to float
for i in range(len(points)):
    points[i] = points[i].split("\t")
    # remove element containing \n
    points[i].pop()
    for j in range(dimensions):
        points[i][j] = float(points[i][j])
        dimension_axis[j].append(points[i][j])

# plot the points
plt.figure()
fig, ax = plt.subplots()

plt.scatter(dimension_axis[0], dimension_axis[1], marker='o', color='b', s=0.5, label='Data Points')


clusters = {
    0: {'centroid': [2.350181, 3.919374], 'count': 831, 'variance': [9707/831 - (1953/831)**2, 16560/831 - (3257/831)**2]},
    1: {'centroid': [5.115044, 2.592920], 'count': 113, 'variance': [3576/113 - (578/113)**2, 1094/113 - (293/113)**2]},
    2: {'centroid': [1.703704, 1.386831], 'count': 243, 'variance': [1574/243 - (414/243)**2, 934/243 - (337/243)**2]},
    8: {'centroid': [-0.834615, 1.671795], 'count': 1560, 'variance': [4057/1560 - (-1302/1560)**2, 9948/1560 - (2608/1560)**2]},
    9: {'centroid': [-0.840909, 5.954545], 'count': 44, 'variance': [73/44 - (-37/44)**2, 1786/44 - (262/44)**2]}
}

# Define a list of colors
colors = ['r', 'g', 'b', 'c', 'm']

# Create plot and add ellipses with different colors
for i, cluster in clusters.items():
    std_deviation_x = math.sqrt(cluster['variance'][0])
    std_deviation_y = math.sqrt(cluster['variance'][1])
    ellipse = patches.Ellipse(cluster['centroid'], std_deviation_x*3, std_deviation_y*3, alpha=0.5, label=f'Cluster {i}', color=colors[i % len(colors)])
    ax.add_patch(ellipse)

# Set aspect to equal to ensure circular shape
# ax.set_aspect('equal', adjustable='box')
padding = 1.0  # adjust as needed
plt.xlim(min(dimension_axis[0]) - padding, max(dimension_axis[0]) + padding)
plt.ylim(min(dimension_axis[1]) - padding, max(dimension_axis[1]) + padding)

plt.title('Data Plot')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
plt.savefig(plotname)
