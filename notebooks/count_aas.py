from collections import Counter
import matplotlib.pyplot as plt
import numpy as np


def count_aas(data):
    with open(data, 'r') as file:
        read_data = file.read()
        read_data_split = read_data.split('\n')

        protein_list = ""
        for line in read_data_split:
            if not line.startswith('>'):
                protein_list += line
        counter = dict(Counter(protein_list))
    return counter


def plot_aa_histogram(counter_as):
    aminoacids = counter_as.keys()
    counted_as = counter_as.values()

    plt.figure(figsize=(10, 8))
    plt.bar(aminoacids, counted_as)
    plt.title('aminoacid distribution')
    plt.xlabel('1-letter-code of aminoacids')
    plt.ylabel('aminoacid count')
    plt.yticks(np.arange(0, 1200000, 50000))
    plt.show()


counter_human = count_aas('uniprot_human.fasta')
print(counter_human)
plot_aa_histogram(counter_human)

counter_animalia = count_aas('uniprot_musmusculus.fasta')
print(counter_animalia)
plot_aa_histogram(counter_animalia)

counter_plantae = count_aas('uniprot_Athaliana.fasta')
print(counter_plantae)
plot_aa_histogram(counter_plantae)

counter_bacteria = count_aas('uniprot_Bacillussubtilis.fasta')
print(counter_bacteria)
plot_aa_histogram(counter_bacteria)

counter_archea = count_aas('uniprot_Methanococcusmaripaludis.fasta')
print(counter_archea)
plot_aa_histogram(counter_archea)