"""Running this script defines and runs the csv standardization function and std save function"""

import pandas as pd

def data_format():
    """Used for standardizing biomarker data with group, timepoint and subject variables in long format.\
    Input csv file path. Input group, subject, timepoint and x_variable names.\
    Group and timepoint variables must contain numeric.\
    Removes rows with more than one field of missing data, standardizes variables, converts time to sorted float."""

    # takes a csv file path and creates data frame. Returns list of column names
    global file_path
    file_path = str(input("Enter csv file path:").strip())
    global file
    file = pd.read_csv(file_path)
    headers = list(file)
    print(headers)

    # inputting column names
    grp_name = 0
    while grp_name not in headers:
        grp_name = str(input("Input group variable name:"))
    tme_name = 0
    while tme_name not in headers:
        tme_name = str(input("Input timepoint variable name:"))
    data_name = 0
    while data_name not in headers:
        data_name = str(input("Input measured data variable name:"))
    subj_name = 0
    while subj_name not in headers:
        subj_name = str(input("Input subject variable name:"))

    # assigns standard names, drops any column that is not group, subject ,time or data
    # drops rows with more than 1 missing value
    file.rename(
        columns={grp_name: "group", subj_name: "subject", tme_name: "time"},
        inplace=True,
    )
    file = file[["group", "time", "subject", str(data_name)]]
    file.dropna(axis=0, thresh=1, inplace=True)

    # convert group, subject, time to string
    file["time"] = file["time"].astype(str)
    file["subject"] = file["subject"].astype(str)
    file["group"] = file["group"].astype(str)

    # extracts only numbers from group and time, removes blanks spaces from subject
    file["time"] = file["time"].str.extract(r"(-?[0-9]+[,./]*[0-9]*)", expand=False)
    file["subject"] = file["subject"].str.replace(" ", "")
    file["group"] = file["group"].str.extract(r"(-?[0-9]+[,./]*[0-9]*)", expand=False)
    
    # changes time to float and sorts by time
    file["time"] = file["time"].astype(float)
    file.sort_values("time")

    # replaces character data in biomarker column with missing value
    file[data_name].replace(r"[a-zA-Z]+", None, inplace=True, regex=True)


def save_std_file():
    """Saves standardized file in the same directory as the input file with _std tag"""
    global file
    global file_path
    sv_file_path = file_path.split(".csv")[0] + "_std.csv"
    file.to_csv(str(sv_file_path), index=False)


# running functions when python file is run from interpreter
if __name__ == "__main__":
    data_format()
    save_std_file()
