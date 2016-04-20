import collections
import numpy as np
import time

from sklearn.decomposition import PCA
from scipy.sparse import csr_matrix


def read_data(file_name):
    num_cols = 10101
    nucleo_dict = {}

    # Initialize dicts for keeping counts of nucleobases
    for c in range(num_cols):
        nucleo_dict[c] = collections.Counter()

    with open(file_name, "r") as f:
        for line in f:
            # Only get nucleobases
            point_entries = line.split(" ")[3:]
            for num, entry in enumerate(point_entries):
                nucleo_dict[num][entry] += 1


    most_common_list = []
    for num, entry_counter in enumerate(nucleo_dict):
        most_common_entry = entry_counter.most_common(1)
        most_common_list.append(most_common_entry[0])

    # the binary matrix that contains real-valued data we will be using
    bin_data_mat = csr_matrix((995, 10101))
    with open(file_name, "r") as f:
        data_num = 0
        for line in f:
            # Only get nucleobases
            point_entries = line.split(" ")[3:]
            for num, entry in enumerate(point_entries):
                if entry != most_common_list[num]:
                    bin_data_mat[data_num, num] = 1.
            data_num += 1

    return nucleo_dict, most_common_list


if __name__ == "__main__":
    nucleo_dict, most_common_list = read_data("p4dataset.txt")

