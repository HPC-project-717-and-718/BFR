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

# clusters = {
#     0: {'centroid': (-0.019072, 2.024158), 'count': 1573, 'variance': [4000/1573 - (-30/1573)**2, 13326/1573 - (3184/1573)**2]},
#     # 1: {'centroid': (4.291720, 0.336306), 'count': 785, 'variance': [69720/785 - (3369/785)**2, 4250/785 - (264/785)**2]},
#     2: {'centroid': (3.318966, 4.909483), 'count': 232, 'variance': [4004/232 - (770/232)**2, 6906/232 - (1139/232)**2]},
#     # 3: {'centroid': (380.500000, 0.000000), 'count': 6, 'variance': [1361517/6 - (2283/6)**2, 73/6 - (0/6)**2]},
#     4: {'centroid': (-2.927326, 2.750000), 'count': 344, 'variance': [4261/344 - (-1007/344)**2, 4938/344 - (946/344)**2]},
#     5: {'centroid': (-2.596774, -0.489247), 'count': 186, 'variance': [2185/186 - (-483/186)**2, 579/186 - (-91/186)**2]},
#     6: {'centroid': (2.365854, -2.506098), 'count': 656, 'variance': [6837/656 - (1552/656)**2, 6312/656 - (-1644/656)**2]},
#     # 7: {'centroid': (8937066.000000, 3.000000), 'count': 2, 'variance': [-2147483648/2 - (17874132/2)**2, 21/2 - (6/2)**2]},
#     # 8: {'centroid': (153329.333333, -0.666667), 'count': 3, 'variance': [-2147483648/3 - (459988/3)**2, 4/3 - (-2/3)**2]},
#     # 7: {'centroid': (6131.750000, 0.750000), 'count': 4, 'variance': [345379178/4 - (24527/4)**2, 34/4 - (3/4)**2]}
# }

# # Define a list of colors
# colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'orange', 'purple', 'lime']

# # Create plot and add ellipses with different colors
# for i, cluster in clusters.items():
#     std_deviation_x = math.sqrt(abs(cluster['variance'][0]))
#     std_deviation_y = math.sqrt(abs(cluster['variance'][1]))
#     ellipse = patches.Ellipse(cluster['centroid'], std_deviation_x*3, std_deviation_y*3, alpha=0.5, label=f'Cluster {i}', color=colors[i % len(colors)])
#     ax.add_patch(ellipse)

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
