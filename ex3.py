# Bioinfotools exercise 3
# Put your code instead of the 'pass' statements
import random
from copy import deepcopy
from Bio import pairwise2
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
from Bio.Align.Applications import MafftCommandline
from io import StringIO
from Bio import AlignIO


class matrix:
    def __init__(self, n_rows, n_cols):
        assert type(n_cols) == int and type(n_rows) == int
        assert n_cols > 0 and n_rows > 0
        self.matrix = [[0 for i in range(n_rows + 1)] for i in range(n_cols + 1)]

    def get(self, i, j):
        return self.matrix[i + 1][j + 1]

    def set(self, x, i, j):
        self.matrix[i + 1][j + 1] = x

    def delete_row(self, i):
        del self.matrix[i + 1]

    def delete_column(self, j):
        for i in range(len(self.matrix) - 1):
            del self.matrix[i + 1][j + 1]

    def get_row(self, i):
        return self.matrix[i + 1]

    def __len__(self):
        return len(self.matrix)


def upgma(D, seq_names_lst):
    seq_names_lst = deepcopy(seq_names_lst)
    weights = {seq: 1 for seq in seq_names_lst}  # As in num of sequences in each cluster.
    while len(seq_names_lst) > 1:
        index1, index2 = find_min(D)
        if index2 < index1:
            index1, index2 = index2, index1  # Even though the matrix is symmetric, work with lower indices.
        recalculate_matrix(D, index1, index2, weights[seq_names_lst[index1]], weights[seq_names_lst[index2]])
        join_clusters(seq_names_lst, index1, index2, weights)
    return seq_names_lst[0]


def find_min(mat):  # mat is a matrix
    min_val = float('inf')
    n = len(mat) - 1
    min_row = -1
    min_column = -1
    for i in range(n):  # Traversing through lower triangle
        for j in range(i):
            if mat.get(i, j) < min_val and mat.get(i, j) != 0:
                min_val = mat.get(i, j)
                min_row = i
                min_column = j
    return min_row, min_column


def join_clusters(seq_names_lst, index1, index2, weights):
    weights[seq_names_lst[index2]] += weights[seq_names_lst[index1]]
    weights.pop(seq_names_lst[index1])  # Rearranging the dictionary
    seq_names_lst[index1] = (seq_names_lst[index1], seq_names_lst[index2])  # Combining two clusters into one.
    weights[seq_names_lst[index1]] = weights.pop(seq_names_lst[index2])
    seq_names_lst.remove(seq_names_lst[index2])


def recalculate_matrix(mat, index1, index2, weight1, weight2):
    for i in range(len(mat) - 1):
        if i != index1 and i != index2:
            x = calculate_average(mat.get(i, index1), mat.get(i, index2), weight1, weight2)
            mat.set(x, i, index1)
            mat.set(x, index1, i)  # Keeping symmetry
    mat.delete_column(index2)
    mat.delete_row(index2)


def calculate_average(num1, num2, weight1, weight2):
    numerator = (num1 * weight1) + (num2 * weight2)
    denominator = weight1 + weight2
    return numerator / denominator


def globalpw_dist(seq_lst):
    n = len(seq_lst)
    dist_mat = matrix(n, n)
    sim_score_mat = matrix(n, n)
    blosum62_matrix = MatrixInfo.blosum62
    max_score = 0
    for j in range(n):  # traverse through lower triangle of matrix - columns
        for i in range(j, n):  # rows
            if i != j:
                curr_score = pairwise2.align.globalds(seq_lst[i], seq_lst[j], blosum62_matrix, -5, -5, score_only=True)
                sim_score_mat.set(curr_score, i, j)
                sim_score_mat.set(curr_score, j, i)
                max_score = max(max_score, curr_score)
    for j in range(n):
        for i in range(j, n):
            if i != j:
                score_cell_value = sim_score_mat.get(i, j)
                dist_mat.set(max_score - score_cell_value + 1, i, j)
                dist_mat.set(max_score - score_cell_value + 1, j, i)
    return dist_mat


def kmer_dist(seq_lst, k=3):
    n = len(seq_lst)
    dist_mat = matrix(n, n)
    kmer_list_per_seq = []
    for i in range(n):
        kmer_list_per_seq.append(calculate_all_k_length_substrings(seq_lst[i], k))
    for i in range(n):
        for j in range(i, n):
            distance = calculate_kmer_dist(kmer_list_per_seq[i], kmer_list_per_seq[j])
            dist_mat.set(distance, i, j)
            dist_mat.set(distance, j, i)  # dist_mat is symmetric
    return dist_mat


def calculate_kmer_dist(k_group_1, k_group_2):
    all_substrings = set(k_group_1.keys())
    seq2_kmers = set(k_group_2.keys())
    all_substrings.update(seq2_kmers)  # Uniting two sets to save memory
    dist = 0
    for word in all_substrings:
        c1 = k_group_1.setdefault(word, 0)
        c2 = k_group_2.setdefault(word, 0)
        dist += (c1 - c2) ** 2
    return dist ** 0.5


def calculate_all_k_length_substrings(seq, k):
    dict_of_substrings = dict()
    for i in range(len(seq) - k + 1):
        temp = seq[i: i+k]
        if temp not in dict_of_substrings:
            dict_of_substrings[temp] = 0
        dict_of_substrings[temp] += 1
    return dict_of_substrings


def eval_dist(seq_lst, msa_aln_path, dist_func=globalpw_dist):
    """
    :param seq_lst: list
        list of n sequences S1, S2, ..., Sn
    :param msa_aln_path: str
        ClustalW FASTA alignment file
    :param dist_func: a distance function name
    :return: dict
        keys are tuples representing each of the splits in the UPGMA tree tree
        values are the counters of the splits in the random permutations
    """
    seq_names, msa_list = split_fasta(msa_aln_path)
    # Second part:
    dist_mat = dist_func(seq_lst)
    tree = upgma(dist_mat, seq_names)
    splits_dict = create_splits_dict(tree)
    for i in range(100):
        columns = choose_random_indices(len(msa_list[0]))
        new_seq_lst = []
        for sequence in msa_list:
            new_seq_lst.append(construct_new_subsequence(sequence, columns))
        new_dist_mat = dist_func(new_seq_lst)
        new_tree = upgma(new_dist_mat, seq_names)
        new_tree_dict = create_splits_dict(new_tree)
        for split in splits_dict:
            x = split in new_tree_dict
            if x:
                splits_dict[split] += 1
    return splits_dict


def choose_random_indices(n):  # This function is not as efficient as possible, but it's mine and O(n) expected
    indices = []
    for i in range(n):
        yes_or_no = random.randint(0, 1)
        if yes_or_no == 1:
            indices.append(i)
    return indices  # Sort this list later to maintain O(n) expected complexity


def construct_new_subsequence(seq, indices):
    sub_list = []
    for index in indices:
        if seq[index] != "-":  # No gaps :)
            sub_list.append(seq[index])
    return "".join(sub_list)


def create_splits_dict(tree):  # tree is a tuple
    list_of_splits = create_splits_tuples(tree)
    return {split: 0 for split in list_of_splits}


def create_splits_tuples(tree):  # tree is a tuple
    list_of_splits = find_splits(tree)
    list_of_splits = clean_singletons_tuplize(list_of_splits)
    return list_of_splits


def find_splits(tree):  # tree is a tuple
    list_of_splits = []
    return find_splits_rec(tree, list_of_splits)


def find_splits_rec(tree, lst):
    if type(tree) == tuple:
        lst.append(sorted([add_tree_items_rec(tree[0], []), add_tree_items_rec(tree[1], [])], key=min))  # Add all leaves of each subtree as tuples.
        lst.sort(key=min)  # Sort the list.
        find_splits_rec(tree[0], lst)  # Apply same logic for left subtree.
        find_splits_rec(tree[1], lst)  # Apply same logic for right subtree.
    return lst


def add_tree_items_rec(tree, splits_lst):
    if type(tree) != tuple:
        splits_lst.append(tree)
        splits_lst.sort(key=min)
    else:
        add_tree_items_rec(tree[0], splits_lst)
        add_tree_items_rec(tree[1], splits_lst)
    return sorted(splits_lst, key=min)


def clean_singletons_tuplize(lst):
    if len(lst) == 1:
        return lst[0]
    for i in range(len(lst)):
        if type(lst[i]) == list:
            lst[i] = clean_singletons_tuplize(lst[i])
    return tuple(lst)


def split_fasta(fasta_path):
    fasta_handle = SeqIO.parse(fasta_path, 'fasta')
    fasta_names = []
    fasta_seqs = []
    for seq_record in fasta_handle:
        fasta_names.append(seq_record.id)
        fasta_seqs.append(seq_record.seq)
    return fasta_names, fasta_seqs


def handle_tree(tree):
    no_apostrophes = str(tree).split("'")
    tree = "".join(no_apostrophes)
    no_spaces = tree.split(" ")
    tree = "".join(no_spaces)
    return tree


def align_sequences(filename, aln_filename):
    mafft_cline = MafftCommandline(MAFFT_EXE_PATH, input=filename)
    stdout, stderr = mafft_cline()
    align = AlignIO.read(StringIO(stdout), "fasta")
    AlignIO.write(align, aln_filename, "fasta")
    return stdout


if __name__ == '__main__':
    MAFFT_EXE_PATH = r"mafft.bat"  # depends on your operating system
    seqs_path = r"sequences.fasta"
    msa_aln_path = r"sequences.aln.fasta"
    # you can write whatever you want here

    # To run section h, please uncomment lines 275-277, and run the program.
    # Please note that running eval_dist with globalpw_dist() might take more than 10 minutes.
    # Remove following lines from comment:

    # seqs = split_fasta(seqs_path)[1]
    # print(eval_dist(seqs, msa_aln_path))
    # print(eval_dist(seqs, msa_aln_path, dist_func=kmer_dist))
    pass
