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
    end_time = time.time()
    print(f"Min seq align value (classic DP): {min_seq_align_val}")
    time_taken = (end_time - start_time)*1000
    return [time_taken, min_seq_align_val]


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
    return divide_and_conquer(generated_strings[0], generated_strings[1])


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
    X: First string from input arguments
    Y: Second string from input arguments
    """
    i = 0
    j = 0
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

    # find the first and second string alignments by reverse engineering the DP optimal result
    max_len = xlen + ylen
    i = xlen
    j = ylen

    # NOTE that these are temporary positions for the current index that we are trying to figure out
    x_pos, y_pos = max_len, max_len
    # create the arrays to store the final aligned strings for x (s1) and y (s2)
    x_aligned, y_aligned = np.zeros(
        max_len + 1, dtype=int), np.zeros(max_len + 1, dtype=int)

    # go through the dp_table (using the 4 cases for which we can reconstruct the string: 1. chars match 2. paying delta cost 3. paying gap cost for s1 4. paying gap cost for s2)
    while i != 0 and j != 0:
        # case 1: chars match
        if x[i - 1] == y[j - 1]:
            x_aligned[x_pos] = ord(x[i - 1])
            y_aligned[y_pos] = ord(y[j - 1])
            i -= 1
            j -= 1
        # case 2: paying the delta cost for mismatch chars
        elif dp_table[i - 1][j - 1] + get_alpha(x[i - 1], y[j - 1]) == dp_table[i][j]:
            x_aligned[x_pos] = ord(x[i - 1])
            y_aligned[y_pos] = ord(y[j - 1])
            i -= 1
            j -= 1
        # case 3: paying the alpha cost for x
        elif dp_table[i][j - 1] + DELTA == dp_table[i][j]:
            x_aligned[x_pos] = ord('_')
            y_aligned[y_pos] = ord(y[j - 1])
            j -= 1
        # case 4: paying the alpha cost for y
        elif dp_table[i - 1][j] + DELTA == dp_table[i][j]:
            x_aligned[x_pos] = ord(x[i - 1])
            y_aligned[y_pos] = ord('_')
            i -= 1
        # decrease the xpos and ypos for each iteration so that we get the current char for x and y
        x_pos -= 1
        y_pos -= 1

    # go through the rest of the string that hasn't been finished yet (since we assumed the answer for each string to be max_len)
    

    while i != 0:
        x_aligned[x_pos] = ord(x[i - 1])
        y_aligned[y_pos] = ord('_')
        i -= 1

    while j != 0:
        x_aligned[x_pos] = ord('_')
        y_aligned[y_pos] = ord(y[j - 1])
        j -= 1


    # remove the extra chars that we didn't use in the x_aligned and y_aligned
    starting_idx = 1  # represents the idx in which the string alignment answer starts at
    i = max_len
    while i >= 1:
        # in this case we have found the start of the junk chars
        if chr(x_aligned[i]) == '_' and chr(y_aligned[i]) == '_':
            starting_idx = i + 1
            break
        i -= 1

    # print out the misaligned genes
    i = starting_idx
    x_final_sequence = ""
    while i <= max_len:
        x_final_sequence += chr(x_aligned[i])
        i += 1
    print(f"s1 aligned: {x_final_sequence}")

    i = starting_idx
    y_final_sequence = ""
    while i <= max_len:
        y_final_sequence += chr(y_aligned[i])
        i += 1
    print(f"s2 aligned: {y_final_sequence}")

    print("S1 Aligned: ", x_final_sequence)
    print("S2 Aligned: ", y_final_sequence)

    # NOTE that we aren't using 0 indexing here because we added a 1-layer padding to the dp_table
    return [x_final_sequence, y_final_sequence, dp_table[xlen][ylen]]



def get_efficient_minimum_penalty(x: str, y: str, is_first: bool):
    """
    X: First string from input arguments
    Y: Second string from input arguments
    is_first: Flag checking if it is running on the first half or last half of string
    """

    xlen = len(x)
    ylen = len(y)

    # initialize dp_table
    dp_table = np.zeros([2, ylen + 1], dtype=int)


    # filling delta values into first row
    for i in range(ylen + 1):
        dp_table[0][i] = DELTA * i



    # Flag: Checking if function is running on first half of sequence
    if is_first:
        for i in range(1, xlen + 1):
            dp_table[1][0] = i * DELTA
            for j in range(1, ylen + 1):
                dp_table[1][j] = min(dp_table[0][j - 1] + get_alpha(x[i - 1], y[j - 1]), 
                                     dp_table[0][j] + DELTA,
                                     dp_table[1][j - 1] + DELTA)
            for j in range(ylen + 1):
                dp_table[0][j] = dp_table[1][j]

    # Flag: Checking if function is running on last half of sequence
    elif not is_first:
        for i in range(1, xlen + 1):
            dp_table[1][0] = i * DELTA
            for j in range(1, ylen + 1):
                dp_table[1][j] = min(dp_table[0][j - 1] + get_alpha(x[xlen - i], y[ylen - j]), 
                                     dp_table[0][j] + DELTA,
                                     dp_table[1][j - 1] + DELTA)

            for j in range(ylen + 1):
                dp_table[0][j] = dp_table[1][j]

    return dp_table[1]


def divide_and_conquer(x: str, y: str):
    """
    X: First string from input arguments
    Y: Second string from input arguments
    """

    # print("X String: " + x)
    # print("Y String: " + y)
    xlen = len(x)
    ylen = len(y)


    xfirst = x[:xlen//2]
    xlast = x[xlen//2:]

    # print("X First Half: " + xfirst)

    # print("X Last Half: " + xlast)

    # divide and conquer algorithm

    # recurrence case once x or y goes down to length 0 or 1
    # print("X Length: ", xlen)
    # print("Y Length: ", ylen)
    
    if xlen < 2 or ylen < 2:
        return get_minimum_penalty(x, y)

    # running get_efficient_minimum_penalty with flag for checking which half of string x
    else:
        first_half = get_efficient_minimum_penalty(xfirst, y, True)
        last_half = get_efficient_minimum_penalty(xlast, y, False)

        # finding ideal split point for string y
        efficient_penalty = [first_half[i] + last_half[ylen - i] for i in range(ylen + 1)]
        print("Efficient String: ", efficient_penalty)
        ysplit = efficient_penalty.index(min(efficient_penalty))
        print("Y Split: ", ysplit)

        print("Y First Half: ", y[:ysplit])
        print("Y Last Half: ", y[ysplit:])
        print("Current First Half X: ", x[:xlen//2])
        print("Current Last Half X: ", x[xlen//2:])
        # recursion back into divide and conquer along half of x and y's ideal split
        get_first = divide_and_conquer(x[:xlen//2], y[:ysplit])
        get_last = divide_and_conquer(x[xlen//2:], y[ysplit:])


        
        # Currently trying to clean out \x00 in the lines below
        

        print(get_first[0])
        print("Get First 0 0: ", get_first[0])
        print("Get First 0 1: ", get_first[1])
        get_first[0].replace('\x00','')
        get_first[1].replace('\x00','')
            

        get_last[0].replace('\x00','')
        get_last[1].replace('\x00','')

        clean_get_first = []
        for i in get_first:
            clean_get_first_new = []
            if (isinstance(i, str)):
                for j in i:
                    clean_get_first_new.append(j.replace('\x00',''))
            clean_get_first.append(clean_get_first_new)


        clean_get_last = []
        for i in get_last:
            clean_get_last_new = []
            if (isinstance(i, str)):
                for j in i:
                    clean_get_last_new.append(j.replace('\x00',''))
            clean_get_last.append(clean_get_last_new)

        # Trying to clean out \x00 until this point

        print("First Half: ", clean_get_first)
        print("Last Half: ", clean_get_last)

    return[get_first[j] + get_last[j] for j in range(3)]

    
def main():
    divide_and_conquer_result = time_wrapper()
    memory_value = process_memory()
    time_value = divide_and_conquer_result[0]
    dandc_value = divide_and_conquer_result[1]
    print(f'Time taken: {time_value} miliseconds')
    print(f'Memory taken: {memory_value} kilobytes')



    # Write to File
    f = open("output.txt",'w')
    f.write(str(dandc_value[0]))
    f.write('\n')
    f.write(str(dandc_value[1]))
    f.write('\n')
    f.write(str(dandc_value[2]))
    f.write('\n')
    f.write("Time Taken: " + str(time_value))
    f.write('\n')
    f.write("Memory Taken: " + str(memory_value))
    f.close()


if __name__ == '__main__':
    main()
