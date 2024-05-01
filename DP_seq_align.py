import numpy as np
import os
import sys
from resource import *
import time
import psutil


# create the hardcoded alpha and delta table
alpha_data = [[0, 110, 48, 94], [110, 0, 118, 48],
              [48, 118, 0, 110], [94, 48, 110, 0]]

ALPHA = np.array(alpha_data)
DELTA = 30


def process_memory():
    process = psutil.Process()
    memory_info = process.memory_info()
    memory_consumed = int(memory_info.rss/1024)
    return memory_consumed


def time_wrapper():
    start_time = time.time()
    min_seq_align_val = start_dp()
    print(f"Min seq align value (classic DP): {min_seq_align_val}")
    end_time = time.time()
    time_taken = (end_time - start_time)*1000
    return time_taken


def generate_string(base, idxs):
    """
    Generate a new string by repeatedly inserting the base string into itself at specified indices.

    Args:
        base (str): The initial string to be modified.
        idxs (list): A list of integers representing the indices at which the base string should be inserted.

    Returns:
        str: The resulting string after all insertions have been performed.

    Examples:
        generate_string("ACTG", [3, 6, 1])
        'ACACTGACTACTGACTGGTGACTACTGACTGG'
        generate_string("TACG", [1, 2, 9)
        'TATTATACGCTATTATACGCGACGCGGACGCG'
    """
    for val in idxs:
        base = base[:int(val) + 1] + base + base[int(val) + 1:]
    return base


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


def start_dp():
    # TODO before submitting double check if we need to ever test with some other directory instead
    # moving into the SampleTestCases
    os.chdir("SampleTestCases")

    # parse command line args
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # parse input file
    with open(input_file, "r") as myfile:
        data = myfile.read().splitlines()

    generated_strings = []
    i = 0
    while i < len(data):
        val = data[i]
        if not val.isnumeric():
            idxs = []
            i += 1
            # parse all the idx values for this base string
            while i < len(data) and data[i].isnumeric():
                idxs.append(data[i])
                i += 1
            # generate the string with this data so far
            generated_strings.append(generate_string(val, idxs))
        else:
            print("Something went wrong during input string generation")
            break

    # start the sequence alignment algorithm
    return get_minimum_penalty(generated_strings[0], generated_strings[1])


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

    # NOTE that we aren't using 0 indexing here because we added a 1-layer padding to the dp_table
    return dp_table[xlen][ylen]


def main():
    print(f'Time taken: {time_wrapper()} miliseconds')
    print(f'Memory taken: {process_memory()} kilobytes')


if __name__ == '__main__':
    main()
