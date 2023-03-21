import pandas as pd
import json 
import collections

# gets the platemap spreadsheet
def get_platemap_data(PLATEMAP_ID, PLATEMAP_SHEET):
    # gets the platemap spreadsheet 
    url = f"https://docs.google.com/spreadsheets/d/{PLATEMAP_ID}/gviz/tq?tqx=out:csv&sheet={PLATEMAP_SHEET}"
    df = pd.read_csv(url, header=None)
    return df