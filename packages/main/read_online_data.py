"""
Be able to read an an spreadsheets online.

author: eablonet
date: #2019

ideas from:
-----------
    https://stackoverflow.com/questions/24063629/loading-a-generic-google-spreadsheet-in-pandas

"""

import pandas as pd
from io import StringIO
import requests


def get_data():
    """Load data online."""
    file = 'https://docs.google.com/spreadsheets/d/'+\
        '1VZ87d4AwZ04ya_wUjRbyaJDenhnoJHOjUEi_X3i-isw'+\
        '/export?format=csv&gid=1642399226'
    act = requests.get(file)
    dataact = act.content.decode('utf-8')

    manip_data = pd.read_csv(
        StringIO(dataact),
        skiprows=[0, 2], decimal=","
    )

    return manip_data
