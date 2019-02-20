import numpy as np
from scipy.spatial import Delaunay
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


def boundary(points):
    points_1 = points.copy()
    points = points[:, :2]
    # hull = ConvexHull(points)
    edges = alpha_shape(points, alpha=0.35, only_outer=True)
    # ax.plot(points[:, 0], points[:, 1], '.')
    for i, j in edges:
        temp1 = points[[i, j], 0]
        temp2 = points[[i, j], 1]
        ax.plot(points[[i, j], 0], points[[i, j], 1])
    plt.show()

    point_set = set()
    footprintPoints = []
    for i, j in edges:
        if i not in point_set:
            footprintPoints.append(points_1[i])
            point_set.add(i)
        if j not in point_set:
            footprintPoints.append(points_1[j])
            point_set.add(j)
    footprintPoints = np.array(footprintPoints)
    Axes3D.plot(ax, xs=footprintPoints[:, 0], ys=footprintPoints[:, 1], zs=footprintPoints[:, 2])
    plt.show()
    return footprintPoints
    # Axes3D.plot(ax, xs=hull.simplices[:, 0], ys=hull.simplices[:, 1], zs=hull.simplices[:, 2])
    # ax.scatter3D(xs=points[:, 0], ys=points[:, 1], zs=points[:, 2])
    """
    Axes3D.plot(xs=points[:,0], ys=points[:,1], zs=points[:,2])
    for simplex in hull.simplices:
        Axes3D.plot(xs=points[simplex, 0], ys=points[simplex, 1], zs=points[simplex, 2])
    plt.plot(points[hull.vertices, 0], points[hull.vertices, 1], points[hull.vertices, 2], 'r--', lw=2)
    plt.plot(points[hull.vertices[0], 0], points[hull.vertices[0], 1], points[hull.vertices[0], 2], 'ro')
    """
    # plt.show()


def shape(points):
    X = points[:, 0]
    Y = points[:, 1]

    some_poly = Polygon(points)

    # Extract the point values that define the perimeter of the polygon
    x, y = some_poly.exterior.coords.xy
    Axes3D.plot(ax, xs=x, ys=y, zs=1)
    plt.show()
    print(x, y)


def alpha_shape(points, alpha, only_outer=True):
    """
    Compute the alpha shape (concave hull) of a set of points.
    :param points: np.array of shape (n,2) points.
    :param alpha: alpha value.
    :param only_outer: boolean value to specify if we keep only the outer border
    or also inner edges.
    :return: set of (i,j) pairs representing edges of the alpha-shape. (i,j) are
    the indices in the points array.
    """
    assert points.shape[0] > 3, "Need at least four points"

    def add_edge(edges, i, j):
        """
        Add a line between the i-th and j-th points,
        if not in the list already
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            assert (j, i) in edges, "Can't go twice over same directed edge right?"
            if only_outer:
                # if both neighboring triangles are in shape, it's not a boundary edge
                edges.remove((j, i))
            return
        edges.add((i, j))

    tri = Delaunay(points)
    edges = set()
    # Loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.vertices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        # Computing radius of triangle circumcircle
        # www.mathalino.com/reviewer/derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle
        a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        s = (a + b + c) / 2.0
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)
        if circum_r < alpha:
            add_edge(edges, ia, ib)
            add_edge(edges, ib, ic)
            add_edge(edges, ic, ia)
    return edges


if __name__ == "__main__":
    test_file = r'C:\Users\Kyler\Documents\GitHub\CSE450Project\capstone\67P_test_scan.csv'
    vertices = np.genfromtxt(test_file, delimiter=",")
    boundary(vertices)
    # shape(vertices)


