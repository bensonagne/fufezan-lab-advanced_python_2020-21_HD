from notebooks.Exercise_2.count_aas import count_aas
from notebooks.Exercise_2.plot_aa_histogram import plot_aa_histogram
from collections import OrderedDict


if __name__ == '__main__':
    counter_human = count_aas('uniprot_human.fasta')
    print(counter_human)
    sorted_as = OrderedDict(sorted(counter_human.items()))
    plot_aa_histogram(sorted_as)
    
    counter_animalia = count_aas('uniprot_musmusculus.fasta')
    print(counter_animalia)
    sorted_as = OrderedDict(sorted(counter_animalia.items()))
    plot_aa_histogram(sorted_as)
    
    counter_plantae = count_aas('uniprot_Athaliana.fasta')
    print(counter_plantae)
    sorted_as = OrderedDict(sorted(counter_plantae.items()))
    plot_aa_histogram(sorted_as)
    
    counter_bacteria = count_aas('uniprot_Bacillussubtilis.fasta')
    print(counter_bacteria)
    sorted_as = OrderedDict(sorted(counter_bacteria.items()))
    plot_aa_histogram(sorted_as)
    
    counter_archea = count_aas('uniprot_Methanococcusmaripaludis.fasta')
    print(counter_archea)
    sorted_as = OrderedDict(sorted(counter_archea.items()))
    plot_aa_histogram(sorted_as)