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
Script to create and populate a repo containing predictions related to the
RSGB 5MHz network for the period May 2009 - February 2015 (inclusive).

Predicted data is stored in the following format

prediction[TX_SITE][RX_SITE][YEAR][MONTH][VOA/P533][RX_PWR] (MUFday voa only)
"""
import json
import os
import re
import subprocess
from tempfile import NamedTemporaryFile

ITURHFPROP_PATH = "/usr/bin/ITURHFProp"
ITURHFPROP_DATA_PATH = "/usr/local/share/iturhfprop/data/"

VOACAP_PATH = "/usr/local/bin/voacapl"
ITSHFBC_PATH = "/home/jwatson/itshfbc"

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


sites = {"GB3RAL":{"lat":51.56, "lng":-1.29, "gain":-0.7},
    "GB3WES":{"lat":54.56, "lng":-2.63, "gain":4.1},
    "GB3ORK":{"lat":59.02, "lng":-3.16, "gain":3.4},
    "G3SET":{"lat":53.39, "lng":-0.57, "gain":4.0},
    "G3WKL":{"lat":52.10, "lng":-0.71, "gain":0.9},
    "G4ZFQ":{"lat":50.73, "lng":-1.29, "gain":4.0},
    "G8IMR":{"lat":50.91, "lng":-1.29, "gain":-31.3},
    "GM4SLV":{"lat":60.29, "lng":-1.43, "gain":-11.2}}

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
                    rx_pwr_list.append(int(line[11:16]))
                m = re.match('^[-\d\s.]+MUFday\s*$', line)
                if m:
                    muf_day_list.append(float(line[11:16]))
        rx_pwr_list.insert(0, rx_pwr_list.pop())
        muf_day_list.insert(0, muf_day_list.pop())

    except Exception as e:
        print(text_in)
        print(e)
    return (rx_pwr_list, muf_day_list)


def get_p553_prediction_df(tx_site, rx_site, year, month):
    #pred_df = pd.DataFrame()
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
    #print(input_file.name)
    #print(output_file.name)
    result_list = []
    try:
        result_buf = []
        with open(output_file.name) as fp:
            for idx, result in enumerate(re.findall('^\d\d,\s+\d\d,\s+[\d.]+,\s+[\d.]+,\s*[-\d.]+\s*$', fp.read(), re.M)):
                result_list.append(float(result.split(',')[4]))
        result_list.insert(0, result_list.pop())
    except Exception as e:
        print(text_in)
        print(e)
    #os.remove(input_file.name)
    #os.remove(output_file.name)
    #pred_df = pd.DataFrame(result_list)
    return result_list

"""Main program starts here"""

prediction_repo = {}

for tx_site in ["GB3RAL", "GB3WES", "GB3ORK"]:
    prediction_repo[tx_site] = {}
    for rx_site in ["G3SET", "G3WKL", "G4ZFQ", "GM4SLV"]:
        prediction_repo[tx_site][rx_site] = {}
        for year in range(2009, 2016):
            prediction_repo[tx_site][rx_site][year] = {}
            for month in range(1,13):
                prediction_repo[tx_site][rx_site][year][month] = {}

                rx_pwr_list = get_p553_prediction_df(tx_site, rx_site, year, month)
                prediction_repo[tx_site][rx_site][year][month]['p533_rx_pwr'] = rx_pwr_list

                rx_pwr_list, muf_day_list = get_voacap_prediction_df(tx_site, rx_site, year, month)
                prediction_repo[tx_site][rx_site][year][month]['voacap_rx_pwr'] = rx_pwr_list
                prediction_repo[tx_site][rx_site][year][month]['voacap_muf_day'] = muf_day_list

with open('predictRepo.json', 'w') as fp:
    json.dump(prediction_repo, fp)
