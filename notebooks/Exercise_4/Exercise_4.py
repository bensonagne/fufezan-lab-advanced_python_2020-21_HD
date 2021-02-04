import requests
from collections import deque
import pandas as pd
import plotly.graph_objects as go


class Protein(object):
    def __init__(self, protein_id):
        self.protein_id = protein_id

    def get_data(self):
        """
        Downloads the fasta-file of a given protein id, extracts the amino-acid sequence and puts it into a list
        :return: list containing the one letter code of the amino acids in the given sequence
        """
        website = 'https://www.uniprot.org/uniprot/' + self.protein_id + '.fasta?fil=reviewed:yes'
        r = requests.get(website)
        fasta_file = self.protein_id + '.fasta'

        with open(fasta_file, 'wb') as file:
            file.write(r.content)
        with open(fasta_file, 'r') as file:
            read_data = file.read()
            read_data_split = read_data.split('\n')

            protein_list = ""
            for line in read_data_split:
                if not line.startswith('>'):
                    protein_list += line
        protein_list = list(protein_list)
        return protein_list

    def map(self, lookup: dict, aa_property: str, window_len=1):
        """

        :param lookup: nested dictionary containing different amino acid properties and the corresponding values
        assigned to the 1-letter amino acid code
        :param aa_property: amino acid property that has to be examined
        :param window_len: length of the sliding window
        :return: list containing the corresponding property of the amino-acid sequence (averaged in the given sliding
    window)
        """
        protein = self.get_data()
        dictionary = lookup[aa_property]
        property_list = [dictionary[aa] for aa in protein]

        window = deque([], window_len)
        average_list = []
        for specific_property in property_list:
            window.append(specific_property)
            average = sum(window) / len(window)
            average_list.append(average)
        return average_list


def get_lookup_dict(csv):
    """
    Makes a nested dictionary out of a given csv-file
    :param csv: given csv-file
    :return: nested dictionary containing different amino acid properties and the corresponding values
    assigned to the 1-letter amino acid code
    """
    aa_properties_dict = pd.read_csv(csv)
    # return pd.Series(data = data['hydropathy index (Kyte-Doolittle method)'].values, index = data['1-letter code'])
    aa_properties_dict = aa_properties_dict.to_dict()

    list_1lettercode = aa_properties_dict['1-letter code'].values()
    lookup_dict = {}
    for pos, aa_property in enumerate(aa_properties_dict.keys()):
        if pos > 2:
            property_dict_list = aa_properties_dict[aa_property].values()
            property_dict = dict(zip(list_1lettercode, property_dict_list))
            lookup_dict[aa_property] = property_dict

    return lookup_dict


def plot_bar_chart(x_values, y_values, title: str):
    """
    Plots a bar chart with the position of the amino acid on the x-axis and the corresponding hydropathy on the y-axis
    and a given title
    :param x_values: list of the position of the amino acids
    :param y_values: list of the hydropathy values a the specific positions
    :param title: Title of the plot
    :return: Bar Chart
    """
    data_fig = [go.Bar(x=x_values, y=y_values, marker_color="white")]
    fig = go.Figure(data=data_fig)
    fig.update_layout(template = 'plotly_dark', title=title, yaxis=dict(title='hydropathy'),
                      xaxis=dict(title='position in the sequence'), font=dict(size=22))
    fig.show()


if __name__ == '__main__':
    p1 = Protein('P32249')
    sequence = p1.get_data()
    lookup_dictionary = get_lookup_dict('amino_acid_properties.csv')
    result_list = p1.map(lookup_dictionary, 'hydropathy index (Kyte-Doolittle method)', 10)
    position_list = list(range(len(result_list)))
    print(result_list)
    plot_bar_chart(position_list, result_list, "Hydropathy plot of P32249 with a sliding window of 10")
