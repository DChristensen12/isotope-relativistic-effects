"""File-name parsing, Google Sheets import, and CSV export. These are I/O-bound, not compute-bound."""

import pandas as pd
import re

try:
    from google.colab import files as _colab_files
except ImportError:
    _colab_files = None  # not running in Colab -- export_dataframe() just saves locally


def file_number(file_name):
    return file_name.split('_')[1].replace('.csv', '') if '_' in file_name else ''


def import_google_spreadsheet_file(URL):
    """Takes the URL of a public Google Sheet ('Anyone with the link' can view) and
    returns it as a DataFrame."""
    # Regular expression to match and capture the necessary part of the URL
    general_url_pattern = r'https://docs\.google\.com/spreadsheets/d/([a-zA-Z0-9-_]+)(/edit#gid=(\d+)|/edit.*)?'

    # Replace function to construct the new URL for CSV export
    get_spreadsheet_ID_and_sheet_ID = lambda m: f'https://docs.google.com/spreadsheets/d/{m.group(1)}/export?' + \
                                               (f'gid={m.group(3)}&' if m.group(3) else '') + 'format=csv'

    New_URL = re.sub(general_url_pattern, get_spreadsheet_ID_and_sheet_ID, URL)

    try:
        df = pd.read_csv(New_URL, low_memory=False, on_bad_lines='skip')
        df_cleaned = df.dropna(how='all')
    except pd.errors.ParserError as e:
        print("Error parsing:", e)
        return pd.DataFrame()

    return df_cleaned


def export_dataframe(data_frame, name):
    """Saves a DataFrame to CSV. Triggers a browser download too, but only when
    actually running in Colab -- the original always called files.download(), which
    only ever worked there anyway."""
    filename = name + '.csv'
    data_frame.to_csv(filename)
    if _colab_files is not None:
        _colab_files.download(filename)
