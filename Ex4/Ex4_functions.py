import pandas as pd
import plotly.graph_objects as go


def get_lookup_dict(csv):
    """
    Makes a nested dictionary out of a given csv-file
    Args:
        csv: csv-file

    Returns:
        nested dictionary containing different amino acid properties and the corresponding values assigned to the
        1-letter amino acid code

    """
    aa_properties_dict = pd.read_csv(csv)
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
    Args:
        x_values: list of the position of the amino acids
        y_values: list of the hydropathy values a the specific positions
        title: Title of the plot

    Returns:
        Bar Chart

    """
    data_fig = [go.Bar(x=x_values, y=y_values, marker_color="white")]
    fig = go.Figure(data=data_fig)
    fig.update_layout(template='plotly_dark', title=title, yaxis=dict(title='hydropathy'),
                      xaxis=dict(title='position in the sequence'), font=dict(size=22))
    fig.show()
