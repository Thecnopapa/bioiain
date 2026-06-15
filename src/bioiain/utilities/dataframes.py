import os, sys, json

import polars as pl



def print_df(df, all_cols=None, all_rows=None):
    if all_cols is None:
        all_cols = True
    if all_rows is None:
        all_rows = False

    if all_cols and all_rows:
        with pl.Config(tbl_cols=df.width ,tbl_rows=len(df)):
            print(df)
    elif all_cols:
        with pl.Config(tbl_cols=df.width):
            print(df)
    elif all_rows:
        with pl.Config(tbl_rows=len(df)):
            print(df)
    return df

