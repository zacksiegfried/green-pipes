import os
import pandas as pd
import pyreadstat as rs
import glob

"""This script should be placed in the same directory as .csv files and run to create .xpt files of the same name"""

dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

# generating list of csv and xpt file paths
csv_paths = glob.glob(dir_path + "/*.csv")
xpt_paths = []
for path in csv_paths:
    xpt_path = os.path.splitext(path)[0] + ".xpt"
    xpt_paths.append(xpt_path)

# creating list of dfs from csv files
df_list = []
for file in csv_paths:
    df = pd.read_csv(file, skipinitialspace = True)
    df_list.append(df)

# removes illegal characters and exports dfs as .xpt files
for i in range(0,len(df_list)):
    frame = df_list[i]
    frame.columns = frame.columns.str.replace(' ', '')
    frame.columns = frame.columns.str.replace("-|\\.|\\/|'|\\[|\\]|\\(|\\)|\\:", '', regex=True)
    export_path = xpt_paths[i]
    rs.write_xport(frame, export_path)
