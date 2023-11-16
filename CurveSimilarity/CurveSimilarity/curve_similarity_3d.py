# Similarity comparison of three-dimensional(3D) curves (consisting of a series of x-, y- and z-coordinates)

import pandas as pd
from curve_similarity import *

def find_procrustes_rotation_matrix(curve1, curve2):
    """
    Find the rotation matrix that makes the two 3D curves most similar to match the best shape.
    Use Singular Value Decomposition (SVD) to determine the best rotation matrix between the two curves.

    Parameters
    ----------
    curve1 : Curve
        Curve 1.
    curve2 : Curve
        Curve 2.

    Returns
    -------
    R : np.array
        The rotation matrix that makes the two curves most similar.

    Raises
    ------
    ValueError
        When the points of the two curves do not agree.

    Singular Value Decomposition (SVD) is a matrix decomposition method, expressed as: A = U * Sigma * V^T,
    where U and V^T are orthogonal matrices representing the left and right singular vectors,
    and Sigma is a diagonal matrix with non-negative singular values on its diagonal.

    In the task of shape alignment, the objective is to find a rotation matrix R that minimizes the Euclidean distance
    between two sets of centralized data. To achieve this, a matrix H is formed as the dot product of the two datasets.
    Performing SVD on H provides U, Sigma, and V^T. The optimal rotation matrix R can be derived as R = V^T * U^T,
    ensuring minimal distance between the datasets after rotation.

    By calculating the H matrix formed by the dot product of the two sets of data and performing the SVD decomposition on it,
    the rotation information can be extracted directly from the decomposition.
    When applied the matrix R formed by this rotation information to one set of data,
    it can be rotated to the direction closest to the other set of data.

    References
    ----------
    [1] Golub, G. H., & Van Loan, C. F. (2013). Matrix computations. JHU press.

    """

    # Both curves must have the same number of points
    if len(curve1) != len(curve2):
        raise ValueError("curve1 and curve2 must have the same length")
    curve1 = np.array(curve1)
    curve2 = np.array(curve2)

    # The curve is already translated in the function procrustes_normalize_curve(curve, estimation_points),
    # so there is no need to center the shape and calculate the translation distance.
    # Calculate the matrix H
    H = np.dot(curve1.T, curve2)
    # Calculate SVD
    U, S, Vt = np.linalg.svd(H)

    # Calculate the rotation matrix
    R = np.dot(Vt.T, U.T)
    # Make sure it's a rotation matrix, not a reflection matrix
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1  # Invert the last vector
        R = np.dot(Vt.T, U.T)
    return R

def apply_rotation_matrix(curve, R):
    """
    Applied matrix transformation to the curve according to the optimal rotation matrix.

    Parameters
    ----------
    curve : Curve
        The curve to be rotated.
    R : np.array
        Rotation matrix.

    Returns
    -------
    Curve
        The curve after matrix transformation.
    """
    return np.dot(curve, R)

def calculate_curve_similarity_3d(curve_standard, curve_estimate, estimation_points, alpha = 1.0):
    """
    Calculate the similarity of two curves, based on the Procrustes normalization and the Fréchet distance.

    Parameters
    ----------
    curve_standard : Curve
        Standard curve.
    curve_estimate : Curve
        The curve to be estimated.
        Different from the 2D algorithm, 3D requires matrix transformation,
        so it is necessary to distinguish the standard curve and the curve to be evaluated.
        Please pay attention to the order of the two curves when passing parameters.
        If only compared the similarity of two curves,
        without the difference between the standard curve and the curve to be evaluated,
        the parameter order of the curves does not matter.
    estimation_points : int
        Number of points of rebalanced curves during Procrustes normalization.
        This number of points determines how many segments the curve should be divided into during the normalization process,
        thus ensuring that both curves have the same number of points before the similarity comparison.
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

    # Use Procrustes normalization to process curves
    normalized_curve_standard = procrustes_normalize_curve(curve_standard, estimation_points)
    normalized_curve_estimate = procrustes_normalize_curve(curve_estimate, estimation_points)

    # Use Singular Value Decomposition (SVD) to determine the best rotation matrix between the two curves
    R = find_procrustes_rotation_matrix(normalized_curve_standard, normalized_curve_estimate)

    # 2D splits the range of [-π, π],
    # and then rotates it at different angles to see which Angle is most similar to the curve.
    # The 3D case is more complicated, and it is difficult to find the optimal rotation angle,
    # but SVD can directly find the rotation matrix that maximizes the similarity between two curves.
    # So here, we directly matrix transform the curve to be evaluated
    # according to the calculated optimal rotation matrix.
    matrix_transformed_curve = apply_rotation_matrix(normalized_curve_estimate, R)

    # Calculate the Fréchet distance between the transformed curve and the standard curve
    frechet_dist = calculate_frechet_distance(normalized_curve_standard, matrix_transformed_curve)

    # Calculate the length of the two curves
    len_curve_standard = calculate_curve_length(normalized_curve_standard)
    len_curve_estimate = calculate_curve_length(normalized_curve_estimate)

    # Compute the relative length difference of the two curves, ranging [0, 1]
    length_diff = abs(len_curve_standard - len_curve_estimate) / (len_curve_standard + len_curve_estimate)
    # Apply weight to the calculated Fréchet distance based on the weight coefficient alpha
    # and the length difference length_diff between the two curves
    weighted_frechet = frechet_dist * (1 + alpha * length_diff)
    # Compute the geometric mean of the lengths of the two curves
    geo_avg_curve_len = np.sqrt(len_curve_standard * len_curve_estimate)
    # Normalize the weighted Fréchet distance with respect to the geometric mean of the curve lengths,
    # ranging [0, 1], to get the final curve similarity
    similarity = max(min(1 - (weighted_frechet / geo_avg_curve_len), 1), 0)

    return similarity

if __name__ == '__main__':
    # Standard flight curve data (that is, the curve for reference),
    # modify the path to the directory where the data file is located.
    df_standard = pd.read_csv('../sample_data/3d/curve_estimate.csv')
    # The flight curve data to be evaluated, modify the path to the directory where the data file is located
    df_estimate = pd.read_csv('../sample_data/3d/curve_standard.csv')

    # Converts the data frame to a list of three-dimensional points (latitude, longitude, height)
    curve_standard = [np.array(point) for point in df_standard[['Longitude', 'Latitude', 'Height']].values.tolist()]
    curve_estimate = [np.array(point) for point in df_estimate[['Longitude', 'Latitude', 'Height']].values.tolist()]

    # Calculate the similarity between two curves, range [0, 1].
	# Set alpha = 5.0. Given the long curve lengths, the Fréchet distance is relatively small compared to the curve length.
    # Therefore, a larger alpha value is set to balance the influence of the curve length during similarity normalization.
    similarity = calculate_curve_similarity_3d(curve_standard, curve_estimate, estimation_points=100, alpha=5.0)
    print("similarity:", similarity)
