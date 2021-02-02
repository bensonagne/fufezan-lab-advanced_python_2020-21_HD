from collections import Counter


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



