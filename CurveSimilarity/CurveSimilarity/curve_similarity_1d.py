# Similarity comparison of one-dimensional(1D) lines (an array of x coordinates)

def levenshtein_distance_recursive(str1, str2):
    """
    Calculate the Levenshtein distance between two 1D lines (an array of x-coordinates).

    Parameters
    ----------
    str1 : list
           List of coordinates for line 1
    str2 : list
           List of coordinates for line 2

    Returns
    -------
    levenshtein_dist : int
            Levenshtein distance
    Levenshtein distance or edit distance:
    The minimum number of single-character edits (insertions, deletions, or substitutions) required to
    transform a string into another string.
    Algorithm Limitations: Using the recursion algorithm to calculate the Levenshtein distance,
    its performance is low, and the time complexity of the algorithm is O(3^n).
    In the future, Dynamic programming can be used to calculate the Levenshtein distance more efficiently
    and avoid redundant calculation.

    References
    ----------
    [1] Levenshtein, V. I. (1966, February). Binary codes capable of correcting deletions, insertions, and reversals.
        In Soviet physics doklady (Vol. 10, No. 8, pp. 707-710).

    """

    # If str1 is empty, the minimum number of operations is to insert all characters of str2
    if len(str1) == 0:
        return len(str2)
    # If str2 is empty, the minimum number of operations is to insert all characters of str1
    elif len(str2) == 0:
        return len(str1)
    # If the two strings are exactly the same, nothing needs to be done and returned 0
    elif str1 == str2:
        return 0
    # Checks whether the last characters of two strings are the same.
    # If the same, d is set to 0 (i.e., no operation is required),
    # otherwise it is set to 1 (i.e., a replacement operation).
    if str1[len(str1) - 1] == str2[len(str2) - 1]:
        d = 0
    else:
        d = 1

    # Use recursion to calculate the Levenshtein distance for the three possible operations and take the minimum:
    # 1. Remove the last character of str1
    # 2. Remove the last character of str2
    # 3. Removes the last character of both
    # (d=0 or d=1, if the last character is different, a substitution operation is required)
    levenshtein_dist = min(levenshtein_distance_recursive(str1, str2[:-1]) + 1,
               levenshtein_distance_recursive(str1[:-1], str2) + 1,
               levenshtein_distance_recursive(str1[:-1], str2[:-1]) + d)
    # Return Levenshtein distance
    return levenshtein_dist

if __name__ == '__main__':
    str1 = [1, 2, 3, 4, 5]
    str2 = [2, 3, 4, 5, 6]

    # Levenshtein distance of sample data str1, str2
    levenshtein_dist = levenshtein_distance_recursive(str1, str2)
    # Normalize Levenshtein distance levenshtein_dist to similarity, range [0, 1]
    similarity_score = 1 - (levenshtein_dist / max(len(str1), len(str2)))

    # levenshtein_dist: 2
    print('levenshtein_dist: ' + str(levenshtein_dist))
    # similarity_score: 0.6
    print('similarity: ' + str(similarity_score))
