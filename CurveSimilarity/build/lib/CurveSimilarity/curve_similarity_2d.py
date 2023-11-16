# Similarity comparison of two-dimensional(2D) curves (consisting of a series of x- and y-coordinates)

import pandas as pd
from curve_similarity import *

def find_procrustes_rotation_angle(curve1, curve2):
    """
    Find the rotation angle that makes the two curves most similar to match the best shape.
    Use cross and dot products to determine the angle difference between two curves.

    Parameters
    ----------
    curve1 : Curve
        Curve 1.
    curve2 : Curve
        Curve 2.

    Returns
    -------
    float
        The rotation angle, in radians, makes the curves most similar.

    Raises
    ------
    ValueError
        When the points of the two curves do not the same.

    References
    ----------
    [1] Gower, J. C. (1975). Generalized procrustes analysis. Psychometrika, 40, 33-51.

    """

    # Both curves must have the same number of points
    if len(curve1) != len(curve2):
        raise ValueError("curve and relativeCurve must have the same length")
    curve1 = np.array(curve1)
    curve2 = np.array(curve2)
    # The sum of the cross products (in 2D) between the corresponding point vectors of the two curves
    numerator = np.sum(curve1[:, 1] * curve2[:, 0] - curve1[:, 0] * curve2[:, 1])
    # The sum of the dot products between the corresponding point vectors of the two curves
    denominator = np.sum(curve1[:, 0] * curve2[:, 0] + curve1[:, 1] * curve2[:, 1])
    # The angle between the accumulated vectors is obtained using the Four-quadrant arctangent function
    return np.arctan2(numerator, denominator)

def rotate_curve(curve: Curve, theta: float) -> Curve:
    """
    Rotate the curve by the given angle.
    The function is implemented using a rotation matrix from linear algebra.

    Parameters
    ----------
    curve : Curve
        The curve to be rotated.
    theta : float
        Rotation Angle (unit: radians): When theta > 0, rotate counterclockwise; When theta < 0, rotate clockwise.
        旋转角度（单位：弧度），theta > 0时，逆时针旋转；theta < 0时，顺时针旋转。

    Returns
    -------
    Curve
        The rotated curve.

    """

    # Define a two-dimensional rotation matrix that rotates points in the plane by angle theta.
    rotation_matrix = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta), np.cos(theta)]
    ])
    # Rotate all the points on the curve
    return [rotation_matrix @ point for point in curve]

def calculate_curve_similarity_2d(curve1, curve2, estimation_points, rotations=101, restrict_rotation_angle = np.pi, alpha = 1.0):
    """
    The similarity of two curves is calculated, based on the Procrustes normalization and the Fréchet distance.

    Parameters
    ----------
    curve1 : Curve
        Curve 1 to be compared.
    curve2 : Curve
        Curve 2 to be compared.
    estimation_points : int
        The number of points used to rebalance the curve during Procrustes normalization.
        It determines how many segments the curve should be divided into during the normalization process,
        thus ensuring that both curves have the same number of points before the similarity comparison.
    rotations : int
        The number of attempts to rotate one curve to find the maximum similarity to another curve.
        The larger the rotations is, the higher the accuracy is, but the less efficient the code runs.
        The rotation Angle of each attempt is the degree of each part
        after dividing [-restrict_rotation_angle, restrict_rotation_angle] into (rotations - 1) parts.
        By default, restrict_rotation_angle is π and rotations are 101,
        which means that dividing the range [-180°, 180°] into 100 parts
        and trying to rotate and compare according to the Angle of each part.
    restrict_rotation_angle : float
        The maximum rotation Angle, which defaults to x (π, 180°),
        in which case the range of rotation angles is [-180°, 180°].
        This parameter can be adjusted to narrow the rotation.
        Increasing this parameter does not make sense, because the rotation angles will be repeated.
    alpha : float
        Weight coefficient of length difference for final similarity evaluation.
        A larger alpha implies that the influence of length difference on the similarity evaluation of the two curves will increase,
        while the impact based on Fréchet distance (which can be understood as shape influence) will decrease.
        alpha = 1.0 means that the influence of length difference and shape difference on similarity is equal.
        Adjust this value up or down as needed.

    Returns
    ----------
    similarity : float
        The similarity between curves in the range [0, 1]

    Raises
    ------
    ValueError
        When restrict_rotation_angle exceeds π.

    Limitations:
    1. Only global similarity is considered but local similarity is not considered.
       For example, both curves A and B are made up of 100 points.
       The last 50 points of A are very similar to the first 50 points of B,
       but the first 50 points of A are very different from the last 50 points of B.
       In this case, the similarity of two curves may be required to be 50%, but the result of this algorithm will be low.
       Dynamic Time Warping (DTW) algorithm can be used to solve this problem in the future.
    2. When comparing similarity, only shape is taken into account without size.
       The Procrustes normalization eliminates the scale differences of the curves,
       which means that two curves with the same shape but very different scales will be considered similar.

    """

    # The maximum value of restrict_rotation_angle is π, beyond which the ValueError is raised
    if abs(restrict_rotation_angle) > np.pi:
        raise ValueError("restrict_rotation_angle cannot be larger than PI")
    # Use Procrustes normalization to process curves
    normalized_curve1 = procrustes_normalize_curve(curve1, estimation_points)
    normalized_curve2 = procrustes_normalize_curve(curve2, estimation_points)

    # Create a list of all possible rotation angles
    thetas_to_check = [0]

    # Find the optimal curve rotation angle obtained by the Procrustes method
    procrustes_theta = find_procrustes_rotation_angle(normalized_curve1, normalized_curve2)
    # Convert the calculated optimal Procrustes rotation angle to the range [-π, π]
    if procrustes_theta > np.pi:
        procrustes_theta = procrustes_theta - 2 * np.pi
    # If the optimal Procrustes rotation angle is not equal to 0 and less than restrict_rotation_angle,
    # add it to the list of rotation angles
    if procrustes_theta != 0 and abs(procrustes_theta) < restrict_rotation_angle:
        thetas_to_check.append(procrustes_theta)

    # The best Procrustes rotation angle is not necessarily the one that minimizes the distance between the two curves.
    # Therefore, more angles need to be examined.
    for i in range(rotations):
        # Distribute the rotation Angle uniformly in rotations - 1 parts
        # over the range of [-restrict_rotation_angle, restrict_rotation_angle]
        theta = -1 * restrict_rotation_angle + (2 * i * restrict_rotation_angle) / (rotations - 1)
        if theta != 0 and theta != np.pi:
            # The resulting thetas_to_check list contains all the rotation angles to check
            thetas_to_check.append(theta)

    # Initialize the Fréchet distance
    min_frechet_dist = float('inf')
    # For each angle, the Fréchet distance between the rotated curve and the standard curve is calculated
    for theta in thetas_to_check:
        # Rotate normalized_curve1
        rotated_curve1 = rotate_curve(normalized_curve1, theta)
        # Calculate the Fréchet distance between the rotated curve and the other curve
        dist = calculate_frechet_distance(rotated_curve1, normalized_curve2)
        if dist < min_frechet_dist:
            # Update the minimum Fréchet distance,
            # and the final result is the Fréchet distance after rotating the curve according to the most accurate rotation angle
            min_frechet_dist = dist

    # Calculate the length of the two curves
    len_curve1 = calculate_curve_length(normalized_curve1)
    len_curve2 = calculate_curve_length(normalized_curve2)

    # Compute the relative length difference of the two curves, ranging [0, 1]
    length_diff = abs(len_curve1 - len_curve2) / (len_curve1 + len_curve2)
    # Apply weight to the calculated Fréchet distance based on the weight coefficient alpha
    # and the length difference length_diff between the two curves
    weighted_frechet = min_frechet_dist * (1 + alpha * length_diff)
    # Compute the geometric mean of the lengths of the two curves
    geo_avg_curve_len = np.sqrt(len_curve1 * len_curve2)

    # Normalize the weighted Fréchet distance with respect to the geometric mean of the curve lengths,
    # ranging [0, 1], to get the final curve similarity
    similarity = max(min(1 - (weighted_frechet / geo_avg_curve_len), 1), 0)

    return similarity

if __name__ == '__main__':
    # Standard flight curve data (that is, the curve for reference),
    # modify the path to the directory where the data file is located.
    df_standard = pd.read_csv('../sample_data/curve_estimate.csv')
    # The flight curve data to be evaluated, modify the path to the directory where the data file is located
    df_estimate = pd.read_csv('../sample_data/curve_standard.csv')

    # Converts the data frame to a list of two-dimensional points (latitude and longitude)
    curve_standard = [np.array(point) for point in df_standard[['Longitude', 'Latitude']].values.tolist()]
    curve_estimate = [np.array(point) for point in df_estimate[['Longitude', 'Latitude']].values.tolist()]

    # Calculate the similarity between two curves, range [0, 1]
    # Set alpha = 5.0. Given the long curve lengths, the Fréchet distance is relatively small compared to the curve length.
    # Therefore, a larger alpha value is set to balance the influence of the curve length during similarity normalization.
    similarity = calculate_curve_similarity_2d(curve_standard, curve_estimate, estimation_points=100, alpha=5.0)
    print("similarity:", similarity)
