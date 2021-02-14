import requests
from collections import deque
import Ex4.Ex4_functions as Ex4_functions


class Protein(object):
    def __init__(self, protein_id):
        """
        Args:
            protein_id: String of the protein id
        """
        self.protein_id = protein_id

    def get_data(self):
        """
        Downloads the fasta-file of a given protein id, extracts the amino-acid sequence and puts it into a list

        Returns:
            list containing the one letter code of the amino acids in the given sequence

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
        returns a list that contains the corresponding property values of the amino acids in the given sequence.

        Args:
            lookup: nested dictionary containing different amino acid properties and the corresponding values assigned to the 1-letter amino acid code
            aa_property: amino acid property that has to be examined
            window_len: length of the sliding window

        Returns:
            list containing the corresponding property of the amino-acid sequence (averaged in the given sliding window)

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


if __name__ == '__main__':
    p1 = Protein('P32249')
    sequence = p1.get_data()
    lookup_dictionary = Ex4_functions.get_lookup_dict('amino_acid_properties.csv')
    print(lookup_dictionary)
    result_list = p1.map(lookup_dictionary, 'hydropathy index (Kyte-Doolittle method)')
    position_list = list(range(len(result_list)))
    print(result_list)
    Ex4_functions.plot_bar_chart(position_list, result_list, "Hydropathy plot of P32249 with a sliding window of 10")
