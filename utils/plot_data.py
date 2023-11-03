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
    points[i] = points[i].split(" ")
    # remove element containing \n
    # points[i].pop()
    for j in range(dimensions):
        points[i][j] = float(points[i][j])
        dimension_axis[j].append(points[i][j])

# plot the points
plt.figure()
fig, ax = plt.subplots()

plt.scatter(dimension_axis[0], dimension_axis[1], marker='o', color='b', s=0.5, label='Data Points')

# Read the cluster data from the file
cluster_data = {}
with open('src/cluster_data.txt', 'r') as f:
    lines = f.readlines()
    for i in range(0, len(lines), 4):  # Adjust for the correct structure of the file
        cluster_index = i  # Extract cluster index
        centroid = [float(val) for val in lines[i].split(":")[-1].strip().split()]
        count = int(lines[i + 1].split(":")[-1].strip().split()[0])
        sum_values = [float(val) for val in lines[i + 2].split(":")[-1].strip().split()]
        sum_squares_values = [float(val) for val in lines[i + 3].split(":")[-1].strip().split()]
        variance = [(sum_squares_values[j] / count) - (sum_values[j] / count) ** 2 for j in range(dimensions)]
        cluster_data[cluster_index] = {'centroid': centroid, 'count': count, 'variance': variance}
# print(cluster_data)

# Define a list of colors
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'orange', 'purple', 'lime']

# Create plot and add ellipses with different colors
for i, cluster in cluster_data.items():
    std_deviation_x = math.sqrt(abs(cluster['variance'][0]))
    std_deviation_y = math.sqrt(abs(cluster['variance'][1]))
    count = cluster['count']
    ellipse = patches.Ellipse(cluster['centroid'], std_deviation_x*3.5, std_deviation_y*3.5, alpha=0.5,
                              label=f'Cluster {i}', color=colors[i % len(colors)])
    ax.add_patch(ellipse)

# Set aspect to equal to ensure circular shape
ax.set_aspect('equal', adjustable='box')

# Adjust the x and y axis limits
padding = 1.0  # adjust as needed
plt.xlim(min(dimension_axis[0]) - padding, max(dimension_axis[0]) + padding)
plt.ylim(min(dimension_axis[1]) - padding, max(dimension_axis[1]) + padding)

plt.title('Data Plot')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
plt.savefig(plotname)
