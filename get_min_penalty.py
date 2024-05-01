import numpy as np

# create the hardcoded alpha and delta table
alpha_data = [[0, 110, 48, 94], [110, 0, 118, 48],
              [48, 118, 0, 110], [94, 48, 110, 0]]

ALPHA = np.array(alpha_data)
DELTA = 30


def get_alpha(x: str, y: str):
    idxs = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return ALPHA[idxs[x]][idxs[y]]


def get_minimum_penalty(x: str, y: str):
    """
    Find Minimum Penalty
    X: X in PDF
    Y: Y in PDF
    alpha: mismatch penalty
    delta: gap penalty
    """

    # Initialize variables
    i = 0
    j = 0

    # Pattern lengths
    xlen = len(x)
    ylen = len(y)

    # Dynamic Programing Table Setup
    dp_table = np.zeros([xlen + 1, ylen + 1], dtype=int)

    dp_table[0: (xlen + 1), 0] = [i * DELTA for i in range(xlen + 1)]
    dp_table[0, 0: (ylen + 1)] = [i * DELTA for i in range(ylen + 1)]

    # Calculate Minimum Penalty
    i = 1
    while i <= xlen:
        j = 1
        while j <= ylen:

            # apply the recurrence relation from classic DP leture (lecture 8)
            dp_table[i][j] = min(dp_table[i - 1][j - 1] + get_alpha(x[i - 1], y[j - 1]),
                                 dp_table[i - 1][j] + DELTA, dp_table[i][j - 1] + DELTA)

            j += 1
        i += 1
