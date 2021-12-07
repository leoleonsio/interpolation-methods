import math
import numpy as np
import scipy.spatial
import startinpy


def distance(p1, p2):
    return math.hypot(p1[0] - p2[0], p1[1] - p2[1])


class Grid:
    def __init__(self, pts, cellsize, nodata_value=-9999):
        self.cellsize = cellsize
        self.pts_3d = pts
        self.nodata = nodata_value
        self.values = []
        x, y, z = zip(*pts)
        pts_2d = list(zip(x, y))
        self.kd = scipy.spatial.KDTree(pts_2d)

        #generate grid for extent
        min_x, min_y, max_x, max_y = min(x), min(y), max(x), max(y)
        row = min_y
        self.nrows, nvalues = 0, 0
        while row <= max_y + cellsize:
            col = min_x
            self.nrows += 1
            grid_row = []
            while col <= max_x + cellsize:
                nvalues += 1

                # centre of the cell makes its coordinates
                grid_row.append([col + cellsize / 2, row + cellsize / 2])
                col += cellsize
            self.values.append(grid_row)
            row += cellsize

        self.ncols = nvalues / self.nrows
        self.min_x = min_x
        self.min_y = min_y

    def interpolate(self, method, *args):
        mapping = {"nn": self.nn,
                   "idw": self.idw,
                   "tin": self.tin,
                   "laplace": self.laplace}

        mapping[method](*args)

    def nn(self):
        for grid_row in self.values:
            for xy in grid_row:
                x, y = xy[0:2]
                d, i = self.kd.query([x, y], k=1)
                height = self.pts_3d[i][2]
                xy.append(height)

    def idw(self, radius1, radius2, power, angle, min_points, max_points):
        bigger_radius = max(radius1, radius2)
        alfa = math.radians(angle)
        for grid_row in self.values:
            for xy in grid_row:
                x, y = xy[0:2]

                # if grid point in samples then take the sample height value
                d, single_nearest = self.kd.query([x, y], k=1)
                if d == 0:
                    xy.append(self.pts_3d[single_nearest][2])
                    continue

                # first filter using the bigger radius because it will use the index, later check if within ellipse
                ids = self.kd.query_ball_point([x, y], bigger_radius)

                points = []
                for i in ids:
                    xp, yp, height = self.pts_3d[i]

                    # source: https://math.stackexchange.com/questions/426150/what-is-the-general-equation-of-the-ellipse-that-is-not-in-the-origin-and-rotate
                    ellipse_equation = ((xp - x) * math.cos(alfa) + (yp - y) * math.sin(alfa)) ** 2 / radius1 ** 2 + \
                                       ((xp - x) * math.sin(alfa) - (yp - y) * math.cos(alfa)) ** 2 / radius2 ** 2

                    if ellipse_equation <= 1:
                        # ellipse contains
                        dist = math.hypot(x - xp, y - yp)
                        weight = dist ** (-power)
                        weighted_height = weight * height
                        points.append((dist, weight, weighted_height))

                if len(points) > max_points > 0:
                    # sort by distance and take max_points closest points
                    points = sorted(points, key=lambda row: row[0])[:max_points]

                if len(points) < min_points or len(points) == 0:
                    xy.append(self.nodata)
                else:
                    _, weight_sum, height_sum = (sum(row) for row in zip(*points))
                    avg_h = height_sum / weight_sum
                    xy.append(avg_h)

    def tin(self):
        def area(tri):
            # triangle edges
            a = distance(tri[0], tri[1])
            b = distance(tri[1], tri[2])
            c = distance(tri[0], tri[2])
            s = (a + b + c) / 2

            #heron's formula
            return math.sqrt(max(0, (s * (s - a) * (s - b) * (s - c))))

        dt = scipy.spatial.Delaunay(self.kd.data)
        for grid_row in self.values:
            for xy in grid_row:
                triangle = dt.simplices[dt.find_simplex(xy)]

                p0, p1, p2 = self.pts_3d[triangle[0]], self.pts_3d[triangle[1]], self.pts_3d[triangle[2]]
                w0, w1, w2 = area((xy, p1, p2)), area((xy, p0, p2)), area((xy, p0, p1))

                avg_h = (w0 * p0[2] + w1 * p1[2] + w2 * p2[2]) / (w0 + w1 + w2)
                xy.append(avg_h)

    def laplace(self):
        dt = startinpy.DT()
        dt.insert(self.pts_3d)
        for row_no, grid_row in enumerate(self.values, start=1):
            for col_no, xy in enumerate(grid_row, start=1):
                result = self.nodata

                # if grid point in samples then take the sample height value
                d, single_nearest = self.kd.query(xy, k=1)
                if d == 0:
                    xy.append(self.pts_3d[single_nearest][2])
                    continue

                new_idx = dt.insert_one_pt(*xy, 0)
                incident_triangles = dt.incident_triangles_to_vertex(new_idx)
                new_points_idx = sorted(set([val for sub in incident_triangles for val in sub]))
                new_points = [dt.get_point(pidx)[0:2] for pidx in new_points_idx]

                # voronoi diagram from the points of incident triangles
                voronoi = scipy.spatial.Voronoi(new_points)
                # indexes of the vertices of the voronoi region for the inserted point
                vertices_idx = voronoi.regions[voronoi.point_region[len(new_points) - 1]]
                region = list([voronoi.vertices[i] for i in vertices_idx])
                region.append(region[0])

                # edges of the voronoi cell
                cell_edges = [(region[i], region[i + 1]) for i in range(len(region) - 1)]

                weight_sum, values_sum = 0, 0
                new_points_idx.remove(new_idx)
                # for each of the edges find the DT point closest to its middle. This will be the corresponding sample point.
                for i, vor_edge in enumerate(cell_edges):
                    min_dist = math.inf
                    weight, height = 0, 0
                    edge_middle = ((vor_edge[0][0] + vor_edge[1][0]) / 2, (vor_edge[0][1] + vor_edge[1][1]) / 2)

                    # minimum distance from the middle of edge to the point means it is the corresponding delaunay edge point
                    for p_idx in new_points_idx:
                        del_point = dt.get_point(p_idx)
                        cur_dist = distance(edge_middle, del_point)

                        if cur_dist < min_dist:
                            # update current values of weight and height
                            min_dist = cur_dist
                            len_voronoi_edge = distance(*vor_edge)
                            len_del_edge = distance(xy, del_point)
                            weight = len_voronoi_edge / len_del_edge
                            height = del_point[2]

                    # use the values for the pair of points with minimum distance
                    weight_sum += weight
                    values_sum += height * weight

                if weight_sum > 0:
                    result = values_sum / weight_sum
                xy.append(result)
                dt.remove(new_idx)

    def save(self, outfile):
        # assign NODATA to pixels out of convex hull
        dt = scipy.spatial.Delaunay(self.kd.data)
        for grid_row in self.values:
            for xyz in grid_row:
                p = xyz[0:2]
                if dt.find_simplex(p) < 0:
                    # point not inside any of the triangles
                    xyz[2] = self.nodata

        with open(outfile, 'w') as f:
            f.write('NCOLS {}\nNROWS {}\nXLLCORNER {}\nYLLCORNER {}\nCELLSIZE {}\nNODATA_VALUE {}\n'
                    .format(self.ncols,
                            self.nrows,
                            self.min_x,
                            self.min_y,
                            self.cellsize,
                            self.nodata))

            # reverse row order for writing the asc file
            for grid_row in self.values[::-1]:
                for _, _, z in grid_row:
                    f.write('{} '.format(z))
                f.write('\n')


def nn_interpolation(list_pts_3d, jparams):
    cellsize = jparams['cellsize']

    g = Grid(list_pts_3d, cellsize)
    g.interpolate('nn')
    g.save(jparams['output-file'])
    print("File written to", jparams['output-file'])


def idw_interpolation(list_pts_3d, jparams):
    cellsize = jparams['cellsize']
    radius1 = jparams['radius1']
    radius2 = jparams['radius2']
    power = jparams['power']
    angle = jparams['angle']
    min_points = jparams['min_points']
    max_points = jparams['max_points']

    g = Grid(list_pts_3d, cellsize)
    g.interpolate('idw', radius1, radius2, power, angle, min_points, max_points)
    g.save(jparams['output-file'])
    print("File written to", jparams['output-file'])


def tin_interpolation(list_pts_3d, jparams):
    cellsize = jparams['cellsize']

    g = Grid(list_pts_3d, cellsize)
    g.interpolate('tin')
    g.save(jparams['output-file'])
    print("File written to", jparams['output-file'])


def laplace_interpolation(list_pts_3d, jparams):
    cellsize = jparams['cellsize']

    g = Grid(list_pts_3d, cellsize)
    g.interpolate('laplace')
    g.save(jparams['output-file'])
    print("File written to", jparams['output-file'])
