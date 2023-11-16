# Common functions

import numpy as np
from scipy.spatial import distance
from curve_types import Point, Curve

def calculate_point_euclidean_distance(point1, point2):
    """
    Calculate the Euclidean Distance between two points to fit points of different dimensions

    Parameters
    ----------
    point1 : Point
             Point 1 coordinate
    point2 : Point
             Point 1 coordinate

    Returns
    -------
    euclidean_distance   double
    The Euclidean Distance between point 1 and point 2

    Raises
    ------
    ValueError
        when the two-point dimensions do not match

    """

    # Convert the input points to a NumPy array
    point1, point2 = np.array(point1), np.array(point2)

    # Check whether the dimensions of the two points match
    if point1.shape != point2.shape:
        raise ValueError("The dimensions of the two points do not match.")

    # Calculate the Euclidean Distance
    euclidean_distance = np.linalg.norm(point1 - point2)
    return euclidean_distance

def point_subtract(v1: Point, v2: Point) -> Point:
    """
    Calculate the difference between two points and get a new vector from v2 to v1
    This function subtracts element-wise between two points,
    i.e., result[0] = v1[0] - v2[0], result[1] = v1[1] - v2[1] ...

    Parameters:
    -----------
    v1 : Point
        Minuend, the end of the vector.
    v2 : Point
        Subtrahend, the start of the vector.

    Returns:
    --------
    Point
    The new point after element-wise subtraction represents a new vector from v2 to v1,
    the size is the distance from v2 to v1.

    """
    return v1 - v2

def calculate_point_euclidean_norm(v: Point) -> float:
    """
    Calculate the Euclidean norm (i.e. length) from the origin to a given point.
    Use the numpy linalg.norm function, which returns the Euclidean norm from the origin to a given point.

    Parameters:
    -----------
    v : Point
        The point to be calculated.

    Returns:
    --------
    float
        Returns the Euclidean norm from the origin to a given point.

    """
    return np.linalg.norm(v)

def extend_point_on_line(p1: Point, p2: Point, dist: float) -> Point:
    """
    Along the lines defined by p1 and p2, look for points at dist distance from p2.
    This function calculates and returns a new point from the direction between p1 and p2,
    starting at p2 and extending out in the direction from p1 to p2.

    Parameters
    ----------
    p1 : Point
        Reference point, which aims to assist in determining the direction of extension
    p2 : Point
        Extended starting point,
        the aim is to extend based on the point in a determined direction to determine the position of the new point.
    dist : float
        The distance to extend, can be a positive or negative number

    Returns
    -------
    Point
        The new point after extending the specified distance in the direction p1 to p2

    """

    # The vector from p1 to p2
    vect = point_subtract(p2, p1)
    # norm: Scale factor
    norm = dist / calculate_point_euclidean_norm(vect)
    # Extend the specified distance based on p2, in the same direction as p1 to p2
    return p2 + norm * vect

def calculate_curve_length(points: Curve) -> float:
    """
    The curve length is calculated by summing the distance between each pair of consecutive points in the list.

    Parameters
    ----------
    points : Curve
        The curve to be calculated.

    Returns
    -------
    float
        The length of the curve

    """

    # zip(list1, list2) function pairs the elements of two (or more) lists according to their index,
    # which results in an ordered pairing of the points on the curve.
    # That is, suppose the points on the curve are a, b, c, d..., zip function gets (a, b), (b, c), (c, d)...
    # Then the Euclidean Distance of each set of paired points is calculated,
    # and summed to obtain the length of the whole curve.
    return sum(calculate_point_euclidean_distance(point1, point2) for point1, point2 in zip(points[:-1], points[1:]))

def resample_curve_uniformly(curve: Curve, num_points: int=100) -> Curve:
    """
    Uniformly resampled curves:
    Given a curve, create a new curve by inserting points such that the curve has num_points (the default is 100).

    Parameters
    ----------
    curve : Curve
        The curve to be resampled.
    num_points : int
        The number of points to target for the new curve, the default is 100.
        This value is used to interpolate the old curve with num_points so that the new curve ends up with num_points.

    Returns
    -------
    Curve
        The new curve with the given number of num_points

    """

    # Calculate the curve length
    curve_len = calculate_curve_length(curve)
    # A curve with num_points is divided into num_points - 1 subsegments,
    # and segment_len is the estimated length of each subsegment.
    segment_len = curve_len / (num_points - 1)
    # Add the first point of the original curve to the list of new curves
    outline_points = [curve[0]]
    # The last point of the original curve
    end_point = curve[-1]
    remaining_curve_points = curve[1:]
    # Since the start and end points are treated separately,
    # we also need to add num_points - 2 points to the curve, so the number of cycles is num_points - 2.
    for _ in range(num_points - 2):
        last_point = outline_points[-1]
        remaining_dist = segment_len
        outline_point_found = False
        while not outline_point_found:
            # Calculate the distance between the last point determined and the next original curve point
            next_point_dist = calculate_point_euclidean_distance(last_point, remaining_curve_points[0])
            # If this distance is less than our expected distance
            if next_point_dist < remaining_dist:
                # Then we reduce the expected distance
                remaining_dist -= next_point_dist
                # Remove the point from the original list of points and set it to the last point,
                # and then the next while loop is passed to rejudge
                last_point = remaining_curve_points.pop(0)
            else:
                # If this distance exceeds the expected distance,
                # extend remaining_dist - next_point_dist distances in the direction last_point points to remaining_curve_points[0].
                # Due to remaining_dist - next_point_dist <= 0,
                # is equivalent to moving in the opposite direction from last_point to remaining_curve_points[0].
                # The new point, the distance from the last determined new curve point is the original estimated length,
                # and the curve shape is as consistent as possible with the original curve
                # because the direction of the vector is taken into account.
                next_point = extend_point_on_line(last_point, remaining_curve_points[0], remaining_dist - next_point_dist)
                # Add next_point to the list of new curves
                outline_points.append(next_point)
                last_point = next_point
                # Enter next for loop
                outline_point_found = True
    # The end point is treated separately
    outline_points.append(end_point)
    return outline_points

def procrustes_normalize_curve(curve, estimation_points):
    """
    Normalized curves by Procrustes Analysis.
    Balance: Insert points into the curve so that the length of each section is as equal as possible,
    which eliminates the effect of the density of curve points.
    Translation: Eliminates the effect of curve position.
    Scaling: Eliminates the effect of curve size.

    Parameters
    ----------
    curve : Curve
        The curve to be normalized.
    estimation_points : int
        Resampling the curve uniformly, the number of final points of the curve
         (the larger this value is, the less efficient the code is, but the more accurate it is).

    Returns
    -------
    Curve
        Normalized curve.

    """

    # Uniformly resampled curve, to ensure that subsequent operations are based on the same number of points
    curve = resample_curve_uniformly(curve, num_points=estimation_points)
    # Calculate the mean of all the points on the curve, and the result is the center of all the points on the curve
    mean = np.mean(curve, axis=0)
    # Translate the curve so that its center is at the origin of the coordinates
    translated_curve = curve - mean
    # Calculate the square root of the mean of the sum of squared distances of each point on the curve from the origin,
    # and this is a measure of how far each point in the curve is from the origin,
    # which helps us understand the size of the curve and normalize it further.
    scale = np.sqrt(np.mean(np.sum(translated_curve**2, axis=1)))
    # Scale the curve so that it has a uniform size
    return translated_curve / scale

def calculate_frechet_distance(P, Q):
    """
    Calculate the Fréchet distance between two curves.
    Use the dynamic programming method, the time complexity is O(nm),
    where n and m are the number of points in the two curves respectively.

    Parameters
    ----------
    P : Curve
        Curve P, one of the curves to be calculated.
    Q: Curve
        Curve Q, another curve to be calculated.
    Both P and Q are sequences of points on a curve.
    Each curve is a numpy array of N x 3, where N is the number of points and 3 represents the x, y, and z coordinates.

    Returns
    ----------
        Fréchet distance.

    Raises
    ------
    ValueError
        Any curve has zero points.

    Fréchet distance: A method of measuring the similarity between two geometric objects,
    especially for comparing the similarity between two curves.
    Specifically, for two curves, Fréchet distance considers all possible paths between two points
    and finds the maximum value of the point-to-point distance on these paths.
    Fréchet distance takes into account the position and order information of curves,
    so it is more nuanced than Hausdorff distance and other curve similarity measures.

    References
    ----------
    [1] Alt, H. and Godau, M., 1995. Computing the Fréchet distance between two polygonal curves.
          International Journal of Computational Geometry & Applications, 5(01n02), pp.75-91.

    """
    n = len(P)
    m = len(Q)

    if n == 0 or m == 0:
        raise ValueError("Input curves are empty.")

    # ca is a cache array to store the distance values that have already been calculated to avoid double calculation.
    ca = -1 * np.ones((n, m))

    def rec(p, q):
        """
        Calculate the Fréchet distance recursively.
        """

        # If the index of the point is less than 0, infinity is returned, indicating unreachable
        if p < 0 or q < 0:
            return float('inf')
        # If the distance between this pair of points has already been calculated, this pair is returned directly
        if ca[p, q] != -1:
            return ca[p, q]
        # The Fréchet distance of the starting point of two curves is the Euclidean Distance between these two points
        if p == 0 and q == 0:
            ca[p, q] = distance.euclidean(P[0], Q[0])
            return ca[p, q]

        # Dynamic programming recursive computation
        elif p == 0:
            ca[p, q] = max(rec(p, q - 1), distance.euclidean(P[p], Q[q]))
        elif q == 0:
            ca[p, q] = max(rec(p - 1, q), distance.euclidean(P[p], Q[q]))
        else:
            ca[p, q] = max(min(rec(p - 1, q), rec(p - 1, q - 1), rec(p, q - 1)), distance.euclidean(P[p], Q[q]))
        return ca[p, q]

    # Returns Fréchet distance of the two curves
    return rec(n - 1, m - 1)
