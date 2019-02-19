import numpy as np
from scipy.spatial import Delaunay
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def boundary(points):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter3D(xs=points[:, 0], ys=points[:, 1], zs=points[:, 2])
    plt.show()
    tri = Delaunay(points)
    print(tri)


if __name__ == "__main__":
    test_file = r"C:\Users\opera\Downloads\67P_test_scan.csv"
    vertices = np.genfromtxt(test_file, delimiter=',')
    boundary(vertices)




