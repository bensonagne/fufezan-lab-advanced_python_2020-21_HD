import sys
import pandas as pd
from math import inf

def msp_to_df(input_file,
              max_seq_len=30,
              min_ce=36, max_ce=40,
              mz_min=135, mz_max=1400,
              keep_sequence_indices=False):
    """
    Function to read spectrum data from .msp file and convert to dataframe.
    Args:
        input_file (str): path to .msp file
        max_seq_len (int): maximum acceptable sequence length
        min_ce (int): minimum collision energy of spectra to be included in df
        max_ce (int): maximum collision energy of spectra to be included in df
        mz_min (int): lower boundary for m/z to be included in df
        mz_max (int): lower boundary for m/z to be included in df
        keep_sequence_indices (boolean): set to True in order to combine both
                                         DataFrames into one

    Returns:
        df (pd.DataFrame or np.array):   spectrum information within defined
                                         parameters [n_spectra, n_features]
        seqs (pd.DataFrame or np.array): sequences (included in df if
                                         keep_sequence_indices is True)
    """
    with open(input_file, 'r') as file:
        def generate_rows():
            skip_spectrum = contains_value = False
            row = {}
            for line in file:
                if line[0] == 'N' and line[1] == 'a':
                    skip_spectrum = contains_value = False
                    row = {}
                    sequence, rest = line[6:].split('/', maxsplit=1)
                    ce = float(rest.rsplit('_', maxsplit=1)[1][:-3])
                    if len(sequence) <= max_seq_len and min_ce <= ce <= max_ce:
                        row['sequence'] = sequence
                    else:
                        skip_spectrum = True
                elif line == '\n':
                    if (not skip_spectrum) and contains_value:
                        yield row
                elif not skip_spectrum and line[0].isdigit():
                    mz, intensity = map(float, line.split('\t')[:2])
                    key = round(mz)
                    if mz_min <= key <= mz_max:
                        contains_value |= intensity > 0
                        row[key] = max(intensity, row.get(key, -inf))
        df = pd.DataFrame.from_records(generate_rows()) \
            .reindex(['sequence'] + list(range(mz_min, mz_max+1)), axis=1) \
            .fillna(0.0)

    if keep_sequence_indices:
        return df.set_index('sequence')
    else:
        return df.drop(columns='sequence'), df.sequence


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python msp_to_df.py <msp_file_name>')
    else:
        file_name = sys.argv[1]
        file_dir = sys.path.dirname(file_name)
        file_prefix = file_name[:-4]

        df, seqs = msp_to_df(file_name)

        df.to_csv(sys.path.join(file_dir, f'{file_prefix}_df.csv'))
        seqs.to_csv(sys.path.join(file_dir, f'{file_prefix}_seqs.csv'))
