import pandas as pd
import os

notebook_dir = os.path.abspath(os.path.dirname(__file__))
package_dir = os.path.abspath(os.path.join(notebook_dir, os.pardir))
#  ^--- use pathlib in your code - yet another module to be covered 
#       maybe in the next advanced course :)
#
pd.set_option("display.max_columns" , 300)

df_high_level = pd.DataFrame(
    data=[
        {'day': 'Monday', 'Topic': 'Intro, recaps and functions'},
        {'day': 'Tuesday', 'Topic': 'Coding philosophy, data flow and some more useful std modules'},
        {'day': 'Wednesday', 'Topic': 'Test driven development, python module, sphinx'},
        {'day': 'Thursday', 'Topic': 'OOP - Object oriented programming'},
        {'day': 'Friday', 'Topic': 'Q&A and code clean up'},
        {'day': '', 'Topic': ''},
        {'day': 'Monday', 'Topic': ''},
        {'day': 'Tuesday', 'Topic': ''},
        {'day': 'Wednesday', 'Topic': ''},
        {'day': 'Thursday', 'Topic': ''},
        {'day': 'Friday', 'Topic': 'Q&A and Tutorium'},


    ]
)

df_details = pd.DataFrame(
    data=[
        {'day': 1, 'Topic': 'Intro'},
        {'day': 1, 'Topic': 'Procedural stuff'},
        {'day': 1, 'Topic': "python basic in 5'"},
        {'day': 1, 'Topic': 'lists and generators'},
        {'day': 1, 'Topic': 'bisect module'},
        {'day': 1, 'Topic': 'Exercise 01 - setup your environment'},
        # ----------------------------
        {'day': 2, 'Topic': 'Functions and dangrous mistakes'},
        {'day': 2, 'Topic': 'Zen of Python and general coding philosophy'},
        {'day': 2, 'Topic': 'csv module'},
        {'day': 2, 'Topic': 'Collections module'},
        {'day': 2, 'Topic': 'Exercise 2'},
        # ----------------------------
        {'day': 3, 'Topic': 'Discussion of Excercises 1 & 2'},
        {'day': 3, 'Topic': 'Basic plotting with plotly'},
        {'day': 3, 'Topic': 'Exercises 3'},
        # -----------------------------
        {'day': 4, 'Topic': 'Discussion of Excercises 3'},
        {'day': 4, 'Topic': "String format"},
        {'day': 4, 'Topic': 'dicts'},
        {'day': 4, 'Topic': 'itertools'},
        {'day': 4, 'Topic': "OOP"},
        # -----------------------------
        # {'day': 3, 'Topic': 'data flow'},
        # {'day': 5, 'Topic': "Basic Python package"},
        {'day': 5, 'Topic': "Test Driven development"},
        {'day': 5, 'Topic': "Auto documentation with Sphinx"},
        # -----------------------------
        {'day': 6, 'Topic': "Pandas 101"},
        {'day': 7, 'Topic': "..."},
        {'day': 8, 'Topic': "..."},
        {'day': 9, 'Topic': "..."},
        {'day': 10, 'Topic': "Tutorium"},
        # -----------------------------
    ]
)

def display_topics(day=1, df=None):
    if df is None:
        df = df_details
    return df[df['day'] == day][['day', 'Topic']].head(20)

def header():
    from IPython.display import display, Markdown

    with open(os.path.join(notebook_dir, 'course_title.md'), 'r') as fh:
        display(Markdown(fh.read()))

