from notebooks.count_aas import count_aas
from notebooks.plot_aa_histogram import plot_aa_histogram

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