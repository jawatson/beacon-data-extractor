#Use this script to create a clean data set prior to running the bextract script.
# This script only needs to be run once.

import pandas as pd

BEACON_CSV = "/home/jwatson/Downloads/selective-beacon-export.csv"
CLEANED_CSV = "beacon_cl.csv"

df = pd.read_csv(BEACON_CSV,
        sep=',',
        header=0,
        usecols=['ts','StnReporting','Year','Month','Day','HourGMT','Minute','GB3RAL','GB3WES','GB3ORK','Noise'])

num_rows = df.shape[0]
print("Number of rows: {:d}".format(df.shape[0]))

df = df[df.ts != "0000-00-00 00:00:00"] # Remove null timestamps

df['ts'] = pd.to_datetime(df['ts'])

df = df[df.ts >= '2009-05-01'] # Calibration started in May 2009

df=df[df.StnReporting.isin(['G3SET','G3WKL','G4ZFQ','G8IMR','GM4SLV'])] # Remove uncalibrated sites


num_clean_rows = df.shape[0]

print("Number of rows (cleaned) {:d} ({:.1f}%):".format(num_clean_rows, (num_clean_rows/num_rows)*100))

df.to_csv(CLEANED_CSV)
