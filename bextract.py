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
This is a short script to automate the comparison of voacap and ITSHFBC
predictions with results obtained from the RSGBs 5MHz beacon project.
"""

import calendar
import datetime
import ephem
import json
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd

matplotlib.style.use('ggplot')
pd.options.mode.chained_assignment = None  # default='warn'

""" CHANGE THESE PARAMETERS TO DEFINE THE CIRCUITS TO EVALUATE """

TX_SITE = 'GB3RAL'
RX_SITE = 'G4ZFQ'
start = (1, 2010) # A (month, year) tuple
#stop = (2, 2015) # (month, year) tuple or None for a single month
#start = (1, 2010) # A (month, year) tuple
stop = None # (month, year) tuple or None for a single month

# External application paths
#BEACON_CSV = "/home/jwatson/Downloads/selective-beacon-export.csv"
BEACON_CSV = "beacon_cl.csv"

DO_PLOTS = True
SAVE_RESIDUALS = False
WINDOW = 60

sites = {"GB3RAL":{"lat":51.56, "lng":-1.29, "gain":-0.7},
    "GB3WES":{"lat":54.56, "lng":-2.63, "gain":4.1},
    "GB3ORK":{"lat":59.02, "lng":-3.16, "gain":3.4},
    "G3SET":{"lat":53.39, "lng":-0.57, "gain":4.0},
    "G3WKL":{"lat":52.10, "lng":-0.71, "gain":0.9},
    "G4ZFQ":{"lat":50.73, "lng":-1.29, "gain":4.0},
    "G8IMR":{"lat":50.91, "lng":-1.29, "gain":-31.3},
    "GM4SLV":{"lat":60.29, "lng":-1.43, "gain":-11.2}}

"""
Sample Data: BeaconID,ts,StnReporting,StnRepQTH,AerialType,Year,Month,Day,HourGMT,Minute,GB3RAL,GB3WES,GB3ORK,Noise,TrackTime,RALHz,WESHz,ORKHz,HeightASL,AerialPolarisation,AerialAlignment

582469,2004-05-03 16:00:00,G4ZFQ,IO90IR,RD,2004,5,3,16,0,70.0599975585938,,,21.41,1.07,873.109985351562,,,8,H,NE/SW
582470,2004-05-03 16:15:00,G4ZFQ,IO90IR,RD,2004,5,3,16,15,69.7200012207031,,,22.68,1.5,873.700012207031,,,8,H,NE/SW
582471,2004-05-03 16:30:00,G4ZFQ,IO90IR,RD,2004,5,3,16,30,73.4899978637695,,,22.6,1.38,874.090026855469,,,8,H,NE/SW
"""

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
"""http://stackoverflow.com/questions/2637293/calculating-dawn-and-sunset-times-using-pyephem"""
def get_greyline_times(site1, site2, month, year):
    sunrises = []
    sunsets = []
    for dom in [1, calendar.monthrange(year, month)[1]]:
        for site in [site1, site2]:
            obs = ephem.Observer()
            obs.date = "{:d}-{:d}-{:d}".format(year, month, dom)
            obs.lat = str(sites[site]['lat'])
            obs.long = str(sites[site]['lng'])
            sunrise = obs.next_rising(ephem.Sun()).datetime()
            sunset = obs.next_setting(ephem.Sun()).datetime()
            sunrises.append(sunrise.hour + sunrise.minute / 60)
            sunsets.append(sunset.hour + sunset.minute / 60)
    """
    For twilight settings;
    obs.horizon = '-6' #-6=civil twilight, -12=nautical, -18=astronomical
    beg_twilight=fred.previous_rising(ephem.Sun(), use_center=True) #Begin civil twilight
    end_twilight=fred.next_setting   (ephem.Sun(), use_center=True) #End civil twilight
    """
    return [(min(sunrises), max(sunrises)), (min(sunsets), max(sunsets))]


"""
Returns the RMS Error.  If a masked array is presented as the first arg, masked
values are ignored.
"""
def get_rmse(l1, l2):
    return np.sqrt(((l1 - l2) ** 2).mean())


def do_analysis(medians_df, voa_pred_df, p533_pred_df):
    # Build a masked array of median values for each hour.  The mask hides missing UTC values.
    medians_ma = ma.masked_values([medians_df['median_pwr'].get(utc, 1.e20) for utc in range(0,24)], 1.e20)
    #print(type(medians_ma))

    p533_corr = ma.corrcoef(medians_ma, np.array(p533_pred_df['rx_pwr']))[0, 1]
    p533_rmse = get_rmse(medians_ma, np.array(p533_pred_df['rx_pwr']))
    voa_corr = ma.corrcoef(medians_ma, np.array(voa_pred_df['rx_pwr']))[0, 1]
    voa_rmse = get_rmse(medians_ma, np.array(voa_pred_df['rx_pwr']))

    voacap_residuals = voa_pred_df['rx_pwr'].subtract(medians_df['median_pwr'])
    p533_residuals = p533_pred_df['rx_pwr'].subtract(medians_df['median_pwr'])

    # OR the mask with a mask for prob muf <= 0.03 (1 day) (True = value is masked.)
    medians_ma.mask = medians_ma.mask | [False if x>0.03 else True for x in voa_pred_df['muf_day'].tolist()]

    muf_day = y = np.array(voa_pred_df['muf_day'].tolist())
    p533_corr_gt_1d = ma.corrcoef(medians_ma, np.array(p533_pred_df['rx_pwr']))[0, 1]
    p533_rmse_gt_1d = get_rmse(medians_ma, np.array(p533_pred_df['rx_pwr']))
    voa_corr_gt_1d = ma.corrcoef(medians_ma, np.array(voa_pred_df['rx_pwr']))[0, 1]
    voa_rmse_gt_1d = get_rmse(medians_ma, np.array(voa_pred_df['rx_pwr']))

    #p533_residuals_gt_1d = p533_pred_df['rx_pwr'].subtract(medians_df['median_pwr'])
    voacap_residuals_gt_1d = np.ma.masked_where(muf_day<=0.03, voacap_residuals).compressed()
    #voacap_residuals_gt_1d = voa_pred_df['rx_pwr'].subtract(medians_df['median_pwr'])
    p533_residuals_gt_1d = np.ma.masked_where(muf_day<=0.03, p533_residuals).compressed()

    voa_residual_mean_gt_1d = np.mean(voacap_residuals_gt_1d)
    voa_residual_sd_gt_1d = np.std(voacap_residuals_gt_1d)
    p533_residual_mean_gt_1d = np.mean(p533_residuals_gt_1d)
    p533_residual_sd_gt_1d = np.std(p533_residuals_gt_1d)

    #print(voacap_residuals)

    return ({"p533_rmse": p533_rmse,
            "p533_corr": p533_corr,
            "voa_rmse": voa_rmse,
            "voa_corr": voa_corr,
            "p533_rmse_gt_1d": p533_rmse_gt_1d,
            "p533_corr_gt_1d": p533_corr_gt_1d,
            "voa_rmse_gt_1d": voa_rmse_gt_1d,
            "voa_corr_gt_1d": voa_corr_gt_1d},
            voacap_residuals,
            p533_residuals,
            voacap_residuals_gt_1d,
            p533_residuals_gt_1d)


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
analysis_df = pd.DataFrame(columns=("month", "p533_rmse", "p533_corr", "voa_rmse", "voa_corr", "voa_residual_mean", "voa_residual_std", "p533_residual_mean", "p533_residual_std"))
voacap_residuals = []
p533_residuals = []
voacap_residuals_gt_1d = []
p533_residuals_gt_1d = []

print('Loading prediction repo')
with open('predictRepo.json') as predict_data_file:
    predict_repo = json.load(predict_data_file)


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
    #medians_df = wdf[(wdf[TX_SITE]>(wdf.Noise+15))].groupby('utc')[[TX_SITE]].median()
    #print (medians_df)
    #medians_df is essentially a list of utc values and medians
    #Method 2a

    medians_df = pd.DataFrame()
    sample_pts = []
    for utc in range(0,24):
        t = datetime.datetime(year,month,15, utc, 0, 0)
        sample_start = (t - datetime.timedelta(minutes=WINDOW/2)).time()
        sample_end = (t + datetime.timedelta(minutes=WINDOW/2)).time()

        sample = wdf.between_time(sample_start, sample_end, include_start=True, include_end=True)[TX_SITE]
        medians_df.set_value(utc, 'median_pwr', sample.median())
        medians_df.set_value(utc, 'median_std', sample.std())
        medians_df.set_value(utc, 'median_sample', sample.shape[0])
        sample_pts.append(sample.tolist())
        #print(sample.tolist())

    p533_pred_df = pd.DataFrame(
        {'rx_pwr': predict_repo[TX_SITE][RX_SITE][str(year)][str(month)]['p533_rx_pwr'],
         'utc': range(0,24),
        })
    voa_pred_df = pd.DataFrame(
        {'rx_pwr': predict_repo[TX_SITE][RX_SITE][str(year)][str(month)]['voacap_rx_pwr'],
         'muf_day': predict_repo[TX_SITE][RX_SITE][str(year)][str(month)]['voacap_muf_day'],
         'utc': range(0,24),
        })

    medians_df = pd.concat([medians_df, p533_pred_df['rx_pwr']], axis=1)
    medians_df.rename(columns = {'rx_pwr':'p533_rx_pwr'}, inplace = True)

    medians_df = pd.concat([medians_df, voa_pred_df[['rx_pwr', 'muf_day']]], axis=1)
    medians_df.rename(columns = {'rx_pwr':'voa_rx_pwr', 'muf_day':'MUFDay'}, inplace = True)

    # Calculate RMS error
    monthly, voacap_monthly_residuals, p533_monthly_residuals, voacap_monthly_residuals_gt_1d, p533_monthly_residuals_gt_1d = do_analysis(medians_df, voa_pred_df, p533_pred_df)

    s_rise, s_set = get_greyline_times(TX_SITE, RX_SITE, month, year)

    #medians_df.to_csv("{:s}_{:s}_{:s}_{:d}_medians.csv".format(TX_SITE, RX_SITE, calendar.month_name[month], year))


    # Do the plot/
    # colours: http://www.flatdesigncolors.com/
    if DO_PLOTS:

        """
        ax = wdf.plot.scatter(x='utc',
                y=TX_SITE,
                xlim=(0,24),
                ylim=(-180, -80),
                color='#54acd2',
                label='Measurement')
        """
        #print(wdf.head)
        plt.figure()
        plt.boxplot(sample_pts, positions=range(0,24))
        ax = plt.gca()

        if not p533_pred_df.empty:
            label="P533 (RMSE={:.2f} r={:.2f} Mean={:.2f} SD={:.2f})".format(float(monthly['p533_rmse']), float(monthly['p533_corr']), np.mean(p533_monthly_residuals), np.std(p533_monthly_residuals))
            ax.plot(p533_pred_df['utc'], p533_pred_df['rx_pwr'], color='#EB6B56', label=label, linewidth=2, marker='o', mec='#EB6B56')

        if not voa_pred_df.empty:
            label= "VOACAP (RMSE={:.2f} r={:.2f} Mean={:.2f} SD={:.2f})".format(float(monthly['voa_rmse']), float(monthly['voa_corr']), np.mean(voacap_monthly_residuals), np.std(voacap_monthly_residuals))
            ax.plot(voa_pred_df['utc'], voa_pred_df['rx_pwr'], color='#00A885', label=label, linewidth=2, marker='o', mec='#00A885')

        ax.get_xaxis().tick_bottom()
        ax.set_ylim([-170, -70])
        ax.get_yaxis().tick_left()
        ax.set_axis_bgcolor('#efefef')
        #ax.set_xticks(np.arange(0,25,1), minor=True)
        #ax.set_xticks(np.arange(0,25,6
        print("Sunrise extends from {:.2f} until {:.2f}".format(s_rise[0], s_rise[1]))
        print("Sunset extends from {:.2f} until {:.2f}".format(s_set[0], s_set[1]))
        ax.axvspan(s_rise[0], s_rise[1], alpha=0.5, color="grey")
        ax.axvspan(s_set[0], s_set[1], alpha=0.5, color="grey")

        ax.set_xlabel("Universal Time (UTC)")
        ax.set_ylabel("Signal Power (dBW)")

        #ax.set_title("{:s} - {:s} {:s} {:d}".format(TX_SITE, RX_SITE, calendar.month_name[month], year))
        ax.legend(loc='upper right')
        #plt.show()
        print("\nVOACAP Monthly Mean: {:.2f} Std.Dev: {:.2f}".format(np.mean(voacap_monthly_residuals), np.std(voacap_monthly_residuals)))
        print("  P533 Monthly Mean: {:.2f} Std.Dev: {:.2f}".format(np.mean(p533_monthly_residuals), np.std(p533_monthly_residuals)))

        plt.savefig("{:s}-{:s}-{:d}-{:d}.png".format(TX_SITE, RX_SITE, year, month))
        plt.close()
    monthly['month'] = datetime.date(year, month, 15)
    analysis_df = analysis_df.append(monthly, ignore_index=True)
    #print(voacap_monthly_residuals)
    voacap_residuals.extend(voacap_monthly_residuals)
    p533_residuals.extend(p533_monthly_residuals)
    voacap_residuals_gt_1d.extend(voacap_monthly_residuals_gt_1d)
    p533_residuals_gt_1d.extend(p533_monthly_residuals_gt_1d)
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
print()
if len(months_list) == 1:
    print("{:s}-{:s} {:s} {:d}".format(TX_SITE, RX_SITE, calendar.month_name[months_list[0][0]], months_list[0][1]))
else:
    print("{:s}-{:s} {:s} {:d} - {:s} {:d} ({:d} months)".format(TX_SITE, RX_SITE, calendar.month_name[start[0]], start[1], calendar.month_name[stop[0]], stop[1], num_rows))
print("Data parsed for {:d} months ({:.2f}%)".format(num_valid_months, 100*(num_valid_months/num_rows)))
print()
print("VOACAP r: {:.2f}  VOACAP RMSE: {:.2f}".format(analysis_df['voa_corr'].mean(), analysis_df['voa_rmse'].mean()))
print("  P533 r: {:.2f}    P533 RMSE: {:.2f}".format(analysis_df['p533_corr'].mean(), analysis_df['p533_rmse'].mean()))
print()
print("VOACAP Mean: {:.2f} Std.Dev: {:.2f}".format(np.mean(voacap_residuals), np.std(voacap_residuals)))
print("  P533 Mean: {:.2f} Std.Dev: {:.2f}".format(np.mean(p533_residuals), np.std(p533_residuals)))
print()
print("MUFday > 1")
print("==========")
print("VOACAP r: {:.2f}  VOACAP RMSE: {:.2f}".format(analysis_df['voa_corr_gt_1d'].mean(), analysis_df['voa_rmse_gt_1d'].mean()))
print("  P533 r: {:.2f}    P533 RMSE: {:.2f}".format(analysis_df['p533_corr_gt_1d'].mean(), analysis_df['p533_rmse_gt_1d'].mean()))
print()
print("VOACAP Mean: {:.2f} Std.Dev: {:.2f}".format(np.mean(voacap_residuals_gt_1d), np.std(voacap_residuals_gt_1d)))
print("  P533 Mean: {:.2f} Std.Dev: {:.2f}".format(np.mean(p533_residuals_gt_1d), np.std(p533_residuals_gt_1d)))
print()
print(len(voacap_residuals))
print(len(voacap_residuals_gt_1d))
if SAVE_RESIDUALS:
    residual_fn = "{:s}_{:s}_voa.rdl".format(TX_SITE, RX_SITE)
    with open(residual_fn, 'w') as residual_file:
        residual_file.writelines(["%s\n" % residual  for residual in voacap_residuals])

    residual_fn = "{:s}_{:s}_p533.rdl".format(TX_SITE, RX_SITE)
    with open(residual_fn, 'w') as residual_file:
        residual_file.writelines(["%s\n" % residual  for residual in p533_residuals])

    residual_fn = "{:s}_{:s}_voa_g1d.rdl".format(TX_SITE, RX_SITE)
    with open(residual_fn, 'w') as residual_file:
        residual_file.writelines(["%s\n" % residual  for residual in voacap_residuals_gt_1d])

    residual_fn = "{:s}_{:s}_p533_g1d.rdl".format(TX_SITE, RX_SITE)
    with open(residual_fn, 'w') as residual_file:
        residual_file.writelines(["%s\n" % residual  for residual in p533_residuals_gt_1d])


analysis_df.to_csv("{:s}_{:s}_summary.csv".format(TX_SITE, RX_SITE))
