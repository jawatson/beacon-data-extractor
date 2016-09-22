#! /usr/bin/env python3
#
# File: bextract.py
#
# Copyright (c) 2016 J.Watson
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.

"""
This is a short script to automate voacap and ITSHFBC predictions
with results obtained from the RSGBs 5MHz beacon project.
"""

import calendar
import datetime
import pandas as pd
import numpy as np
import numpy.ma as ma
import math
import matplotlib.pyplot as plt
import matplotlib
import os
import re
import io
import subprocess
from tempfile import NamedTemporaryFile

import sunrise

matplotlib.style.use('ggplot')
pd.options.mode.chained_assignment = None  # default='warn'

""" CHANGE THESE PARAMETERS TO DEFINE THE CIRCUITS TO EVALUATE """


TX_SITE = 'GB3RAL'
RX_SITE = 'G4ZFQ'
#start = (7, 2005) # A (month, year) tuple
#stop = (6, 2012) # (month, year) tuple or None for a single month
start = (1, 2010) # A (month, year) tuple
stop = None # (month, year) tuple or None for a single month

# External application paths
#BEACON_CSV = "/home/jwatson/Downloads/selective-beacon-export.csv"
BEACON_CSV = "beacon_cl.csv"
ITURHFPROP_PATH = "/usr/bin/ITURHFProp"
ITURHFPROP_DATA_PATH = "/home/jwatson/github/proppy/flask/data/"
VOACAP_PATH = "/usr/local/bin/voacapl"
ITSHFBC_PATH = "/home/jwatson/itshfbc"
DO_PLOTS = True


sites = {"GB3RAL":{"lat":51.56, "lng":-1.29, "gain":-0.7},
    "GB3WES":{"lat":54.56, "lng":-2.63, "gain":4.1},
    "GB3ORK":{"lat":59.02, "lng":-3.16, "gain":3.4},
    "G3SET":{"lat":53.39, "lng":-0.57, "gain":4.0},
    "G3WKL":{"lat":52.10, "lng":-0.71, "gain":0.9},
    "G4ZFQ":{"lat":50.73, "lng":-1.29, "gain":4.0},
    "G8IMR":{"lat":50.91, "lng":-1.29, "gain":-31.3},
    "GM4SLV":{"lat":60.29, "lng":-1.43, "gain":-11.2}}

ssn_repo = {"sources": ["http://sidc.oma.be/silso/INFO/snmstotcsv.php", "http://sidc.oma.be/silso/FORECASTS/prediSC.txt"],
    "retrieved": 1472203198.672534,
    "ssn": {
        "2009": {"2": 2.7, "3": 2.9, "10": 10.9, "12": 12.7, "11": 11.7, "5": 3.5, "9": 9.5, "4": 3.3, "7": 5.5, "1": 2.5, "6": 4.1, "8": 7.4},
        "2010": {"2": 16.1, "3": 18.5, "10": 34.5, "12": 42.5, "11": 39.1, "5": 23.1, "9": 29.5, "4": 20.8, "7": 25.2, "1": 14.0, "6": 24.6, "8": 26.4},
        "2011": {"2": 48.8, "3": 53.8, "10": 87.4, "12": 92.5, "11": 89.4, "5": 69.3, "9":86.6, "4": 61.1, "7": 83.6, "1": 45.7, "6": 77.2, "8": 86.3},
        "2012": {"2": 98.2, "3": 98.3, "10": 85.8, "12": 88.1, "11": 87.7, "5": 90.9, "9": 85.3, "4": 95.1, "7": 84.5, "1": 95.5, "6": 86.6, "8": 85.1},
        "2013": {"2": 86.1, "3": 84.4, "10": 107.0, "12": 107.6, "11": 106.9, "5": 87.0, "9": 104.7, "4": 84.3, "7": 94.6, "1": 86.8, "6": 90.9, "8": 99.0},
        "2014": {"2": 110.5,"3": 114.3, "10":97.3, "12":92.2, "11":94.7, "5":115.0, "9":101.9, "4":116.4, "7":112.6, "1":109.3, "6":114.1, "8":108.3},
        "2015": {"2": 86.1, "3": 82.2, "10": 64.3, "12": 57.8, "11": 61.3, "5": 76.1, "9": 65.9, "4": 78.9, "7": 68.3, "1": 89.3, "6": 72.1, "8": 66.4},
        "2016": {"2": 51.8, "3": 48.0, "10": 28.1, "12": 25.5, "11": 26.8, "5": 41.2, "9": 30.2, "4": 44.4, "7": 35.1, "1": 54.5, "6": 38.0, "8": 32.7},
        "2017": {"4": 21.2, "3": 22.1, "6": 19.8, "7": 19.1, "2": 22.9, "5": 20.6, "1": 24.1}}}


voa_ssn_repo = {"sources": ["ftp://ftp.ngdc.noaa.gov/STP/space-weather/solar-data/solar-indices/sunspot-numbers/predicted/table_international-sunspot-numbers_monthly-predicted.txt",],
    "retrieved": 1472203198.672534,
    "ssn": {
        "2009": {"1":1.8, "2":1.9,   "3":2.0,  "4":2.2,  "5":2.3,  "6":2.7,  "7":3.6,  "8":4.8,  "9":6.1,  "10":7.1,  "11":7.6,  "12":8.3},
        "2010": {"1":9.3, "2":10.6,  "3":12.3, "4":14.0, "5":15.5, "6":16.4, "7":16.8, "8":17.4, "9":19.6, "10":23.2, "11":26.5, "12":28.8},
        "2011": {"1":30.9, "2":33.4, "3":36.9, "4":41.8, "5":47.6, "6":53.2, "7":57.2, "8":59.0, "9":59.5, "10":59.9, "11":61.1, "12":63.4},
        "2012": {"1":65.5, "2":66.9, "3":66.8, "4":64.6, "5":61.7, "6":58.9, "7":57.8, "8":58.2, "9":58.1, "10":58.6, "11":59.7, "12":59.6},
        "2013": {"1":58.7, "2":58.4, "3":57.6, "4":57.9, "5":59.9, "6":62.6, "7":65.5, "8":68.9, "9":73.0, "10":74.9, "11":75.3, "12":75.9},
        "2014": {"1":77.3, "2":78.4, "3":80.8, "4":81.9, "5":80.5, "6":79.7, "7":78.6, "8":75.6, "9":70.8, "10":67.3, "11":65.4, "12":63.7},
        "2015": {"1":61.9, "2":60.5, "3":59.7, "4":59.1, "5":57.8, "6":55.5, "7":53.1, "8":51.2, "9":49.7, "10":48.6, "11":47.4, "12":46.1}}}

"""
Sample Data: BeaconID,ts,StnReporting,StnRepQTH,AerialType,Year,Month,Day,HourGMT,Minute,GB3RAL,GB3WES,GB3ORK,Noise,TrackTime,RALHz,WESHz,ORKHz,HeightASL,AerialPolarisation,AerialAlignment

582469,2004-05-03 16:00:00,G4ZFQ,IO90IR,RD,2004,5,3,16,0,70.0599975585938,,,21.41,1.07,873.109985351562,,,8,H,NE/SW
582470,2004-05-03 16:15:00,G4ZFQ,IO90IR,RD,2004,5,3,16,15,69.7200012207031,,,22.68,1.5,873.700012207031,,,8,H,NE/SW
582471,2004-05-03 16:30:00,G4ZFQ,IO90IR,RD,2004,5,3,16,30,73.4899978637695,,,22.6,1.38,874.090026855469,,,8,H,NE/SW
"""

def get_voacap_prediction_df(tx_site, rx_site, year, month):
    tx_lat = float(sites[tx_site]['lat'])
    tx_lat_str = "{:.2f}{:s}".format(abs(tx_lat), 'N' if tx_lat>=0 else 'S')
    tx_lng = float(sites[tx_site]['lng'])
    tx_lng_str = "{:.2f}{:s}".format(abs(tx_lng), 'E' if tx_lng>=0 else 'W')
    tx_gain = float(sites[tx_site]['gain'])
    rx_lat = float(sites[rx_site]['lat'])
    rx_lat_str = "{:.2f}{:s}".format(abs(rx_lat), 'N' if rx_lat>=0 else 'S')
    rx_lng = float(sites[rx_site]['lng'])
    rx_lng_str = "{:.2f}{:s}".format(abs(rx_lng), 'E' if rx_lng>=0 else 'W')
    rx_gain = float(sites[rx_site]['gain'])
    ssn = voa_ssn_repo['ssn'][str(year)][str(month)]
    buf = []

    buf.append('LINEMAX      55')
    buf.append('COEFFS    CCIR')
    buf.append('TIME          1   24    1    1')
    buf.append('LABEL     {:<20s}{:<20s}'.format(tx_site, rx_site))
    buf.append('CIRCUIT{:>9s}{:>10s}{:>10s}{:>10s}  S     0'.format(tx_lat_str, tx_lng_str, rx_lat_str, rx_lng_str))
    buf.append('SYSTEM       1. 150. 3.00  90. 38.0 3.00 0.10')
    buf.append('FPROB      1.00 1.00 1.00 0.00')
    buf.append('ANTENNA       1    1   02   30{:>10.3f}[samples/sample.00    ]  0.0    0.0090'.format(tx_gain))
    buf.append('ANTENNA       2    2   02   30     0.000[samples/sample.00    ]  0.0{:>10.4f}'.format(rx_gain))
    buf.append('FREQUENCY  5.29')
    buf.append('METHOD       30    0')
    buf.append('MONTH      {:d}{:>5.2f}'.format(year, month))
    buf.append('SUNSPOT  {:>6.1f}'.format(float(ssn)))
    buf.append('EXECUTE')
    buf.append('QUIT')
    text_in = "{:s}\n".format('\n'.join(buf))
    input_fn = "nvis.dat"
    output_fn = "nvis.out"

    with open(os.path.join(*[ITSHFBC_PATH, "run", input_fn]), "w") as input_file:
        input_file.write(text_in)

    FNULL = open(os.devnull, 'w')

    return_code = subprocess.call([VOACAP_PATH,
        ITSHFBC_PATH,
        input_fn,
        output_fn],
        stdout=FNULL,
        stderr=subprocess.STDOUT)
    result_list = []
    rx_pwr_list = []
    muf_day_list = []
    try:
        with open(os.path.join(*[ITSHFBC_PATH, "run", output_fn])) as fp:
            for line in fp:
                m = re.match('^[-\d\s]+S DBW\s*$', line)
                if m:
                    rx_pwr_list.append(line)
                m = re.match('^[-\d\s.]+MUFday\s*$', line)
                if m:
                    muf_day_list.append(line)
        for idx in range(0, 24, 1):
            rx_pwr = int(rx_pwr_list[idx][11:16])
            muf_day = float(muf_day_list[idx][11:16])
            result_list.append({"utc":idx+1, "rx_pwr":rx_pwr, "muf_day":muf_day})
        midnight = result_list[-1]
        result_list.insert(0, {"utc":0, "rx_pwr":midnight['rx_pwr'], "muf_day":midnight['muf_day']})
    except Exception as e:
        print(text_in)
        print(e)
    pred_df = pd.DataFrame(result_list)
    return pred_df


def get_p553_prediction_df(tx_site, rx_site, year, month):
    pred_df = pd.DataFrame()
    tx_lat = float(sites[tx_site]['lat'])
    tx_lng = float(sites[tx_site]['lng'])
    tx_gain = float(sites[tx_site]['gain'])
    rx_lat = float(sites[rx_site]['lat'])
    rx_lng = float(sites[rx_site]['lng'])
    rx_gain = float(sites[rx_site]['gain'])
    ssn = ssn_repo['ssn'][str(year)][str(month)]
    buf = []
    buf.append('Path.L_tx.lat {:.2f}'.format(tx_lat))
    buf.append('Path.L_tx.lng {:.2f}'.format(tx_lng))
    buf.append('TXAntFilePath "{:s}"'.format('ISOTROPIC'))
    buf.append('TXGOS {:.2f}'.format(tx_gain))
    buf.append('RXAntFilePath "{:s}"'.format('ISOTROPIC'))
    buf.append('RXGOS {:.2f}'.format(rx_gain))
    buf.append('Path.year {:d}'.format(year))
    buf.append('Path.month  {:d}'.format(month))
    buf.append('Path.hour 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24')
    buf.append('Path.SSN {:.2f}'.format(ssn))
    buf.append('Path.frequency 5.29')
    buf.append('Path.txpower {:.2f}'.format(-20.458)) # 9W
    buf.append('Path.BW {:.1f}'.format(500))
    buf.append('Path.SNRr {:.1f}'.format(0))
    buf.append('Path.SNRXXp 90')
    buf.append('Path.ManMadeNoise "{:s}"'.format('RURAL'))
    buf.append('Path.Modulation "ANALOG"')
    buf.append('Path.SorL "SHORTPATH"')
    buf.append('RptFileFormat "RPT_OPMUF | RPT_PR"')
    buf.append('LL.lat {:.2f}'.format(rx_lat))
    buf.append('LL.lng {:.2f}'.format(rx_lng))
    buf.append('LR.lat {:.2f}'.format(rx_lat))
    buf.append('LR.lng {:.2f}'.format(rx_lng))
    buf.append('UL.lat {:.2f}'.format(rx_lat))
    buf.append('UL.lng {:.2f}'.format(rx_lng))
    buf.append('UR.lat {:.2f}'.format(rx_lat))
    buf.append('UR.lng {:.2f}'.format(rx_lng))
    buf.append('DataFilePath "{:s}"'.format(ITURHFPROP_DATA_PATH))
    input_file = NamedTemporaryFile(mode='w+t', prefix="proppy_", suffix='.in', delete=False)
    text_in = "{:s}\n".format('\n'.join(buf))
    input_file.write(text_in)
    input_file.close()

    FNULL = open(os.devnull, 'w')
    output_file = NamedTemporaryFile(prefix="proppy_", suffix='.out', delete=False)
    return_code = subprocess.call([ITURHFPROP_PATH,
        input_file.name,
        output_file.name],
        stdout=FNULL,
        stderr=subprocess.STDOUT)
    print(input_file.name)
    #print(output_file.name)
    result_list = []
    try:
        result_buf = []
        with open(output_file.name) as fp:
            for idx, result in enumerate(re.findall('^\d\d,\s+\d\d,\s+[\d.]+,\s+[\d.]+,\s*[-\d.]+\s*$', fp.read(), re.M)):
                result_list.append({"utc":idx+1, "rx_pwr":float(result.split(',')[4])})
        midnight = result_list[-1]
        result_list.insert(0, {"utc":0, "rx_pwr":midnight['rx_pwr']})
    except Exception as e:
        print(text_in)
        print(e)
    #os.remove(input_file.name)
    os.remove(output_file.name)
    pred_df = pd.DataFrame(result_list)
    return pred_df


def get_utc(row):
    return row['HourGMT'] + (row['Minute'] / 60.0)


"""
Calibration levels supplied by Dr. M.Walden.
NOTE: Calibration levels are valid from May 2009 only
"""
def calibrate_rx_level(rx_db, rx_site):
    offset = 0
    if rx_site == "G3WKL":
        offset = 168.1
    elif rx_site == "G3SET":
        offset = 180.00
    elif rx_site == "G4ZFQ":
        offset = 177.9
    elif rx_site == "GM4SLV":
        offset = 183
    elif rx_site == "G8IMR":
        offset = 179
    else:
        raise KeyError("Calibration not found for site {:s}".format(rx_site))
    return rx_db - offset


"""
Returns a pair of tuples containing sun rise/set activity
for the link.  Each tuple contains the earliest / latest
occurance for the event chosen from the min/max of the
event times for the start and end of the month at each
end of the link.
"""
def get_greyline_times(site1, site2, month, year):
    sunrises = []
    sunsets = []

    for dom in [1, calendar.monthrange(year, month)[1]]:
        for site in [site1, site2]:
            s = sunrise.sun(lat=float(sites[site]['lat']), long=float(sites[site]['lng']))
            dt = datetime.datetime(year, month, dom)
            st = s.sunrise(when=dt)
            sunrises.append(st.hour + st.minute / 60)
            st = s.sunset(when=dt)
            sunsets.append(st.hour + st.minute / 60)
    return [(min(sunrises), max(sunrises)), (min(sunsets), max(sunsets))]


def get_correlation(l1, l2):
    return np.corrcoef(l1, l2)[0, 1]


def get_rmse(l1, l2):
    return np.sqrt(((l1 - l2) ** 2).mean())


def do_analysis(medians_df, sample_sizes, voa_pred_df, p533_pred_df):
    # Build a masked array of median values for each hour.  The mask hides missing UTC values.
    medians_ma = ma.masked_values([medians_df[TX_SITE].get(utc, 1.e20) for utc in np.arange(0,24,1)], 1.e20)

    p533_corr = get_correlation(medians_ma, np.array(p533_pred_df['rx_pwr'].tolist()[:-1]))
    p533_rmse = get_rmse(medians_ma, np.array(p533_pred_df['rx_pwr'].tolist()[:-1]))
    voa_corr = get_correlation(medians_ma, np.array(voa_pred_df['rx_pwr'].tolist()[:-1]))
    voa_rmse = get_rmse(medians_ma, np.array(voa_pred_df['rx_pwr'].tolist()[:-1]))

    # OR the mask with a mask for prob muf <= 0.03 (1 day) (True = value is masked.)
    medians_ma.mask = medians_ma.mask | [False if x>0.03 else True for x in voa_pred_df['muf_day'].tolist()[:-1]]

    # Calculate RMS and P between predicted and measured values
    p533_corr_gt_1d = get_correlation(medians_ma, np.array(p533_pred_df['rx_pwr'].tolist()[:-1]))
    p533_rmse_gt_1d = get_rmse(medians_ma, np.array(p533_pred_df['rx_pwr'].tolist()[:-1]))
    voa_corr_gt_1d = get_correlation(medians_ma, np.array(voa_pred_df['rx_pwr'].tolist()[:-1]))
    voa_rmse_gt_1d = get_rmse(medians_ma, np.array(voa_pred_df['rx_pwr'].tolist()[:-1]))

    return {"p533_rmse": p533_rmse, "p533_corr": p533_corr,
            "voa_rmse": voa_rmse, "voa_corr": voa_corr,
            "p533_rmse_gt_1d": p533_rmse_gt_1d, "p533_corr_gt_1d": p533_corr_gt_1d,
            "voa_rmse_gt_1d": voa_rmse_gt_1d, "voa_corr_gt_1d": voa_corr_gt_1d}


def get_null_entry(year, month):
    return {"month":datetime.date(year, month, 15),
        "p533_rmse": np.nan, "p533_corr": np.nan,
        "voa_rmse": np.nan, "voa_corr": np.nan}


"""
Generates a list of month/year tuples from specified start stop dates.
e.g. start=(2,2015) stop=(6,2016) expands to;
[(2,2015), (3,2015), (4,2015) ... (4,2016),(5,2016),(6,2016)]
"""
def get_months_list(start, stop):
    ml = [start]
    if stop is not None:
        tmp = start
        while (tmp[0] < stop[0]) or (tmp[1] < stop[1]):
            m, y = tmp
            m = m+1
            if m > 12:
                m = 1
                y = y+1
            tmp = (m,y)
            ml.append(tmp)
    return ml


# START OF MAIN APPLICATION

#Dataframe to store the results
analysis_df = pd.DataFrame(columns=("month", "p533_rmse", "p533_corr", "voa_rmse", "voa_corr"))

print("Loading RSGB Dataframe...")

# NOTE: The df has a non-unique index
df = pd.read_csv(BEACON_CSV,
        sep=',',
        header=0,
        index_col='ts',
        parse_dates = True,
        usecols=['ts','StnReporting','Year','Month','Day','HourGMT','Minute','GB3RAL','GB3WES','GB3ORK','Noise'])

months_list = get_months_list(start,stop)

for month, year in months_list:
    # Grab the relevent subset of the df and perform the prediction
    print("{:s} - {:s} {:s} {:d}".format(TX_SITE, RX_SITE, calendar.month_name[month], year))

    wdf = df[(df.Year == year) & (df.Month == month) & (df.StnReporting == RX_SITE)]

    if wdf.empty:
        print("No data found for {:s}-{:s}-{:d}-{:d}".format(TX_SITE, RX_SITE, year, month))
        analysis_df = analysis_df.append(get_null_entry(year, month), ignore_index=True)
        continue

    # Test for rx. offline
    if len(wdf.index) <= 1500:
        print("WARNING: Size of dataframe is {:d}".format(len(wdf.index)))

    wdf[TX_SITE] = wdf[TX_SITE].apply(calibrate_rx_level, args=(RX_SITE,))
    wdf['Noise'] = wdf['Noise'].apply(calibrate_rx_level, args=(RX_SITE,))

    wdf['utc'] = wdf.apply(lambda row: get_utc(row), axis=1)
    wdf = wdf[(wdf[TX_SITE]>(wdf.Noise+15))] # Discard values where SNR <= 15

    # Calculate median for values where s/n > 15
    # A List of the sample size at each utc value
    samples_size_list = wdf[(wdf[TX_SITE]>(wdf.Noise+15))].groupby('utc')[[TX_SITE]].count()[TX_SITE].tolist()
    if np.sum(samples_size_list)<1000:
        print("WARNING: Tx. Appears to be offline - skipping this one")
        analysis_df = analysis_df.append(get_null_entry(year, month), ignore_index=True)
        continue
    #Method 1
    medians_df = wdf[(wdf[TX_SITE]>(wdf.Noise+15))].groupby('utc')[[TX_SITE]].median()
    #print (medians_df)
    #medians_df is essentially a list of utc values and medians
    #Method 2
    for window in [0, 30, 60]:
        for utc in range(0,24):
            t = datetime.datetime(year,month,15, utc, 0, 0)
            sample_start = (t - datetime.timedelta(minutes=window/2)).time()
            sample_end = (t + datetime.timedelta(minutes=window/2)).time()

            sample = wdf.between_time(sample_start, sample_end, include_start=True, include_end=True)[TX_SITE]
            #print(sample)
            #print("{:d} Sample size = {:d}.  Median = {:.6f}".format(utc, sample.shape[0], sample.median()))
            medians_df.set_value(utc, 'median_w{:d}'.format(window), sample.median())
            medians_df.set_value(utc, 'median_w{:d}_std'.format(window), sample.std())
            medians_df.set_value(utc, 'median_w{:d}_ss'.format(window), sample.shape[0])
    #print (medians_df)
    medians_df.to_csv("{:s}_{:s}_medians.csv".format(TX_SITE, RX_SITE))

    p533_pred_df = get_p553_prediction_df(TX_SITE, RX_SITE, year, month)
    voa_pred_df = get_voacap_prediction_df(TX_SITE, RX_SITE, year, month)

    # Calculate RMS error
    monthly = do_analysis(medians_df, samples_size_list, voa_pred_df, p533_pred_df)

    s_rise, s_set = get_greyline_times(TX_SITE, RX_SITE, month, year)

    # Do the plot
    # colours: http://www.flatdesigncolors.com/
    if DO_PLOTS:
        ax = wdf.plot.scatter(x='utc',
                y=TX_SITE,
                xlim=(0,24),
                ylim=(-180, -80),
                color='#54acd2',
                label='Measurement')

        medians_df.index.name = 'utc'
        #print(medians_df.head())
        medians_df.reset_index(inplace=True)
        medians_df.plot.scatter(x='utc', y=TX_SITE, color='#2C82C9', label='Median', s=75, ax=ax)

        """ Some temporary plots follow"""
        label="Median W0"
        medians_df.plot.line(x='utc', y='median_w0', label=label, linewidth=2, marker='o', ax=ax)

        label="Median W30"
        medians_df.plot.line(x='utc', y='median_w30', label=label, linewidth=2, marker='o', ax=ax)

        label="Median W60"
        medians_df.plot.line(x='utc', y='median_w60', label=label, linewidth=2, marker='o', ax=ax)

        """end of the temporary plots"""

        if not p533_pred_df.empty:
            label="Predicted (P533) (RMSE={:.2f} r={:.2f})".format(float(monthly['p533_rmse']), float(monthly['p533_corr']))
            p533_pred_df.plot.line(x='utc', y='rx_pwr', color='#EB6B56', label=label, linewidth=2, marker='o', mec='#EB6B56', ax=ax)

        if not voa_pred_df.empty:
            label= "Predicted (VOA) (RMSE={:.2f} r={:.2f})".format(float(monthly['voa_rmse']), float(monthly['voa_corr']))
            voa_pred_df.plot.line(x='utc', y='rx_pwr', color='#00A885', label=label, linewidth=2, marker='o', mec='#00A885', ax=ax)

        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ax.set_axis_bgcolor('#efefef')
        ax.set_xticks(np.arange(0,25,1), minor=True)
        ax.set_xticks(np.arange(0,25,6))
        ax.axvspan(s_rise[0], s_rise[1], alpha=0.5, color="grey")
        ax.axvspan(s_set[0], s_set[1], alpha=0.5, color="grey")
        ax.set_xlabel("Universal Time (UTC)")
        ax.set_ylabel("Signal Power (dBW)")
        ax.set_title("{:s} - {:s} {:s} {:d}".format(TX_SITE, RX_SITE, calendar.month_name[month], year))
        #plt.show()
        plt.savefig("{:s}-{:s}-{:d}-{:d}.png".format(TX_SITE, RX_SITE, year, month))
        plt.close()
    monthly['month'] = datetime.date(year, month, 15)
    analysis_df = analysis_df.append(monthly, ignore_index=True)
"""
print(analysis_df)
col_names = ["p533_corr_gt_1d", "p533_rmse_gt_1d", "voa_corr_gt_1d", "voa_rmse_gt_1d"]
for col in col_names:
    print("{:s} (Mean): {:.2f}".format(col, analysis_df[col].mean()))
"""
analysis_df.voa_rmse_gt_1d = analysis_df.voa_rmse_gt_1d.astype(float)
analysis_df.p533_rmse_gt_1d = analysis_df.p533_rmse_gt_1d.astype(float)


#Create a summary plot
if stop:
    ax = analysis_df.plot.line(x='month',
            y='voa_rmse_gt_1d',
            ylim=(0,40),
            figsize=(10,3),
            color='#00A885',
            label='VOACAP',
            linewidth=2)

    analysis_df.plot.line(x='month',
            y='p533_rmse_gt_1d',
            color='#EB6B56',
            label="P533",
            linewidth=2,
            ax=ax)

    ax.set_title("{:s}-{:s} {:s} {:d} - {:s} {:d}".format(TX_SITE, RX_SITE, calendar.month_name[start[0]], start[1], calendar.month_name[stop[0]], stop[1]))

    plt.savefig("{:s}_{:s}_summary.png".format(TX_SITE, RX_SITE))

#End of main loop
num_rows = len(analysis_df.index)
num_valid_months = analysis_df['p533_corr_gt_1d'].count()
print(analysis_df)
print()
if len(months_list) == 1:
    print("{:s}-{:s} {:s} {:d}".format(TX_SITE, RX_SITE, calendar.month_name[months_list[0][0]], months_list[0][1]))
else:
    print("{:s}-{:s} {:s} {:d} - {:s} {:d} ({:d} months)".format(TX_SITE, RX_SITE, calendar.month_name[start[0]], start[1], calendar.month_name[stop[0]], stop[1], num_rows))
print("Data parsed for {:d} months ({:.2f}%)".format(num_valid_months, 100*(num_valid_months/num_rows)))
print()
print("  P533 r: {:.2f}    P533 RMSE: {:.2f}".format(analysis_df['p533_corr'].mean(), analysis_df['p533_rmse'].mean()))
print("VOACAP r: {:.2f}  VOACAP RMSE: {:.2f}".format(analysis_df['voa_corr'].mean(), analysis_df['voa_rmse'].mean()))
print()
print("MUFday > 1")
print("==========")
print("  P533 r: {:.2f}    P533 RMSE: {:.2f}".format(analysis_df['p533_corr_gt_1d'].mean(), analysis_df['p533_rmse_gt_1d'].mean()))
print("VOACAP r: {:.2f}  VOACAP RMSE: {:.2f}".format(analysis_df['voa_corr_gt_1d'].mean(), analysis_df['voa_rmse_gt_1d'].mean()))

analysis_df.to_csv("{:s}_{:s}_summary.csv".format(TX_SITE, RX_SITE))
