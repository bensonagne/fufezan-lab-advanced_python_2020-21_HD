import pandas as pd
from collections import deque
import plotly.graph_objects as go


def extract_sequence_from_fasta(fasta):
    """
    Extracts the amino-acid sequence from the fasta file into a list
    Args:
        fasta: fasta file

    Returns:
        list containing the one letter code of the amino acids in the given sequence

    """
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
    """
    Makes a dictionary from a csv files with the column '1-letter code' as keys and the column
    'hydropathy index (Kyte-Doolittle method)' as values.
    Args:
        csv: csv file

    Returns:
        dictionary with the 1-letter code as keys and the corresponding hydropathy as values

    """
    data = pd.read_csv(csv)
    # return pd.Series(data = data['hydropathy index (Kyte-Doolittle method)'].values, index = data['1-letter code'])
    data = data.to_dict()
    list_1lettercode = data['1-letter code'].values()
    list_hydropathy = data['hydropathy index (Kyte-Doolittle method)'].values()
    dictionary = dict(zip(list_1lettercode, list_hydropathy))
    return dictionary


def hydropathy_value_list_of_sequence(sequence, mapping_dict, window_len=None):
    """
    returns a list that contains the corresponding hydropathy values of the amino acids in the given sequence.
    Args:
        sequence: list of the amino-acid sequence that should be compared
        mapping_dict: dictionary with the 1-letter code as keys and the corresponding hydropathy as values
        window_len: length of the sliding window

    Returns:
        list containing the corresponding hydropathy of the amino-acid sequence (averaged in the given sliding window)

    """
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


def plot_bar_chart(x_values, y_values, title: str):
    """
    Plots a bar chart with the position of the amino acid on the x-axis and the corresponding hydropathy on the y-axis
    and a given title
    Args:
        x_values: list of the position of the amino acids
        y_values: list of the hydropathy values a the specific positions
        title: Title of the plot

    Returns:
        Bar Chart

    """
    data_fig = [go.Bar(x=x_values, y=y_values, marker_color="white")]
    fig = go.Figure(data=data_fig)
    fig.update_layout(template='plotly_dark', title=title,
                      yaxis=dict(title='hydropathy'), xaxis=dict(title='position in the sequence'))
    fig.show()


if __name__ == '__main__':
    first = extract_sequence_from_fasta('gpcr.fasta')
    second = dictionary_from_csv('amino_acid_properties.csv')
    result = hydropathy_value_list_of_sequence(first, second, 20)
    position_list = list(range(len(result)))
    print(result)
    plot_bar_chart(position_list, result, "Hydropathy plot with a sliding window of 20")
