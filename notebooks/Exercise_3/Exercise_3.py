import pandas as pd
from collections import deque
import plotly.graph_objects as go


def extract_sequence_from_fasta(fasta):
    with open(fasta, 'r') as file:
        # return [protein_list for line in file.readlines() if not line.startswith('>') for protein_list in
        # line.strip()]
        read_data = file.read()
        read_data_split = read_data.split('\n')

        protein_list = ""
        for line in read_data_split:
            if not line.startswith('>'):
                protein_list += line
    protein_list = list(protein_list)
    return protein_list


def dictionary_from_csv(csv):
    data = pd.read_csv(csv)
    # return pd.Series(data = data['hydropathy index (Kyte-Doolittle method)'].values, index = data['1-letter code'])
    data = data.to_dict()
    list_1lettercode = data['1-letter code'].values()
    list_hydropathy = data['hydropathy index (Kyte-Doolittle method)'].values()
    dictionary = dict(zip(list_1lettercode, list_hydropathy))
    return dictionary


def hydropathy_value_list_of_sequence(sequence, mapping_dict, window_len=None):
    hydropathy_list = [mapping_dict[aa] for aa in sequence]

    if window_len is not None:
        window = deque([], window_len)
        average_list = []
        for hydropathy in hydropathy_list:
            window.append(hydropathy)
            average = sum(window) / len(window)
            average_list.append(average)
        hydropathy_list = average_list
    return hydropathy_list


# plot bar charts
first = extract_sequence_from_fasta('gpcr.fasta')
second = dictionary_from_csv('amino_acid_properties.csv')
hydropathy_liste = hydropathy_value_list_of_sequence(first, second)
print(hydropathy_liste)
position_list = list(range(361))

data_fig = [go.Bar(x=position_list, y=hydropathy_liste)]
fig = go.Figure(data=data_fig)
fig.update_layout(title="Hydropathy plot without a sliding window",
                  yaxis=dict(title='hydropathy'), xaxis=dict(title='position in the sequence'))
fig.show()
