import os


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


# TODO before submitting double check if we need to ever test with some other directory instead
# moving into the SampleTestCases
os.chdir("SampleTestCases")

# take file name as input
filename = input("Enter the input file name: ")

# parse input file
with open(filename, "r") as myfile:
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
