import collections
import matplotlib.pyplot as plt
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

    individual_to_pop = []
    individual_to_gender = []

    with open(file_name, "r") as f:
        count = 0
        for line in f:
            # Only get nucleobases
            point_entries = line.split(" ")[1:]
            gender = point_entries[0]
            population = point_entries[1]
            individual_to_pop.append(population)
            individual_to_gender.append(gender)
            genome = point_entries[2:]
            for num, entry in enumerate(genome):
                nucleo_dict[num][entry] += 1
            count += 1

    # Get distinct entries as list
    distinct = list(set(individual_to_pop))

    # Convert to indexed array
    indexed_pop = [distinct.index(x) for x in individual_to_pop]


    # Get distinct entries as list
    distinct_gender = list(set(individual_to_gender))

    # Convert to indexed array
    indexed_gender = [distinct_gender.index(x) for x in individual_to_gender]


    most_common_list = []
    for _, entry_counter in nucleo_dict.iteritems():
        most_common_entry = entry_counter.most_common(1)
        most_common_list.append(most_common_entry[0])

    # the binary matrix that contains real-valued data we will be using
    bin_data_mat = np.zeros((995, 10101))
    #print most_common_list
    with open(file_name, "r") as f:
        data_num = 0
        for line in f:
            #print "New example processed"
            # Only get nucleobases
            point_entries = line.split(" ")[3:]
            for num, entry in enumerate(point_entries):
                if entry != most_common_list[num][0]:
                    bin_data_mat[data_num, num] = 1.
            data_num += 1

    return nucleo_dict, most_common_list, bin_data_mat, indexed_pop, distinct, indexed_gender, distinct_gender


def run_PCA(X, n_comp=2):
    # Run PCA on X with desired num principle components
    pca = PCA(n_components=n_comp)
    pca.fit(X)

    return pca.transform(X)


def plot_data(components, pc,  pop_indices, pop_names, x_axis, y_axis, title):
    plt.scatter(components[:, pc[0]], components[:, pc[1]], c=pop_indices)
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.title(title)
    plt.show()


def identify_strongest_components(arr, n_comp):
    sorted_arr = np.argsort(arr)
    return sorted_arr[-n_comp:]


if __name__ == "__main__":
    nucleo_dict, most_common_list, bin_data_mat, pop_indices, pop_names, indexed_gender, distinct_gender = read_data("p4dataset.txt")
    print pop_names
    #pca_2d = run_PCA(bin_data_mat, n_comp=2)
    pca_3d = run_PCA(bin_data_mat, n_comp=3)
    #plot_data(pca_2d, pop_indices, pop_names, "V1", "V2", "First Two Principal Components")
    #plot_data(pca_3d, [0, 2], indexed_gender, pop_names, "V1", "V3", "First and Third Principal Components")
    print "Strongest: ", identify_strongest_components(pca_3d[:, 2], 3)

