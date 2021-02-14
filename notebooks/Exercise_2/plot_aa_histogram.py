import matplotlib.pyplot as plt
import numpy as np


def plot_aa_histogram(counter_as):
    aminoacids = counter_as.keys()
    counted_as = counter_as.values()

    plt.figure(figsize=(10, 8))
    plt.bar(aminoacids, counted_as)
    plt.title('aminoacid distribution')
    plt.xlabel('1-letter-code of aminoacids')
    plt.ylabel('aminoacid count')
    plt.show()
