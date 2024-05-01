import numpy as np


def get_minimum_penalty(x: str, y: str, alpha: int, delta: int):
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

    dp_table[0: (xlen + 1), 0] = [i * delta for i in range(xlen + 1)]
    dp_table[0, 0: (ylen + 1)] = [i * delta for i in range(ylen + 1)]

    # Calculate Minimum Penalty
    i = 1
    while i <= xlen:
        j = 1
        while j <= ylen:

			# case 1 (when the chars match up)
            if x[i - 1] == y[j - 1]:
                dp_table[i][j] = dp_table[i - 1][j - 1]

			# case 2 (when chars don't match or misaligned)
            else:
				# apply the recurrence relation from classic DP leture (lecture 8)
                dp_table[i][j] = min(dp_table[i - 1][j - 1] + alpha,
                                     dp_table[i - 1][j] + delta, dp_table[i][j - 1] + delta)

            j += 1
        i += 1