import sys
import math
from Bio import SeqIO
import matplotlib.pyplot as plt


# populates a dictionary where keys are all possible valid codons and values
# are 0
def get_codons():
    initialized_dict = {}
    possible_letters = ['A', 'C', 'G', 'T']
    for i in range(len(possible_letters)):
        for j in range(len(possible_letters)):
            for k in range(len(possible_letters)):
                possible_codon = possible_letters[i] + possible_letters[j] + possible_letters[k]
                initialized_dict[possible_codon] = 0

    return initialized_dict

# returns true if given string is a valid codon else false
def check_codon_validity(codon):
    if (len(codon) != 3):
        return False
    for i in range(len(codon)):
        if codon[i] != "A" and codon[i] != "C" and codon[i] != "G" and codon[i] != "T":
            return False
        
    return True

# returns a populated dictionary with codon counts given a gene sequence
def count_codon(seq):
    codon_dict = get_codons()
    for i in range(0, len(seq), 3):
        cur_codon = seq[i:i+3]
        #print(cur_codon)
        if check_codon_validity(cur_codon):
            codon_dict[cur_codon] += 1
    
    return codon_dict

# returns distance value given the spike and window codon counts
def get_distance(spike_dict, window_dict):
    ns = 0
    nw = 0
    for x in spike_dict:
        ns += spike_dict[x]
        nw += window_dict[x] 
    
    d = 0
    for x in window_dict:
        spike_codon_prop = spike_dict[x]/ns
        window_codon_prop = window_dict[x]/nw
        diff = spike_codon_prop - window_codon_prop
        d += (diff * diff)
    
    return d

 
#iteratively computes codon counts for 2000 bp windows in full genome
def compute_window_distances(genome_seq, spike_dict):
    window_distances = []

    for i in range(len(genome_seq)-2000+1):
        cur_window = genome_seq[i:i+2000]
        cur_window_dict = count_codon(cur_window)
        window_distances.append((i, get_distance(spike_dict, cur_window_dict)))

    return window_distances

#prints window distances to a file given distances list and filename           
def print_distances_to_file(distances, filename):
    with open(filename, 'w') as file:
        for i in range(len(distances)):
            (pos, value) = distances[i]
            file.write(str(pos) + "\t".expandtabs(4) + str(value))
            file.write("\n")

#prints codon counts to a file given a dictionary and filename
def print_dict_to_file(dict, filename):
    with open(filename, 'w') as file:
        for x in dict:
            file.write(x + "\t".expandtabs(4) + str(dict[x]))
            file.write("\n")

#special print for the purpose of producing x and y values in 2 different files
def special_print(list):
    with open("xvalues.txt", 'w') as file:
        for i in range(len(list)):
            file.write(str(i))
            file.write("\n")
    with open("yvalues.txt", 'w') as file:
        for i in range(len(list)):
            (pos, val) = list[i]
            file.write(str(val))
            file.write("\n")

if __name__ == "__main__":

    #reading in spike gene
    record = SeqIO.read("spike.fasta", "fasta")
    gene = str(record.seq)

    #reading in full genome
    record2 = SeqIO.read("fullgenome.fasta", "fasta")
    full_genome = str(record2.seq)

    spike_codon_counts = count_codon(gene)
    
    print_dict_to_file(spike_codon_counts, "spike_codon_count.txt")

    window_distances = compute_window_distances(full_genome, spike_codon_counts)
    print_distances_to_file(window_distances, "window_distance_table.txt")

    #for producing graph purposes
    special_print(window_distances)
