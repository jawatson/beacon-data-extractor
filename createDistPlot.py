"""
This script produces a graph showing the distribution of values from a list
of values contained withing the nominated files.

The default values below bin the errors into 6dB chunks, centred about zero.
"""
import numpy as np
import matplotlib.pyplot as plt
import sys

def rmse(rdls):
    return np.sqrt(np.mean([i ** 2 for i in rdls]))
    #return np.sqrt(sum([i ** 2 for i in residuals]) / len(residuals))


# Constants
BIN_WIDTH = 6
CLIP_LOWER = -50
CLIP_UPPER = 50

values = []
"""
if len(sys.argv) != 2:
    print("USAGE:")
    print("python3 createDifferenceCSV.py difference_table.csv")
    sys.exit(1)
"""
residuals = []
for file_number in range(2, len(sys.argv)+1):
    with open(sys.argv[file_number-1]) as residual_file:
        content = [float(x.strip('\n')) for x in residual_file.readlines()]
    residuals = residuals + content

#Calculate the mean and std.dev before clipping...
count = len(residuals)
mean = np.mean(residuals)
std_dev = np.std(residuals)
rmse = rmse(residuals)

residuals = np.clip(residuals, CLIP_LOWER, CLIP_UPPER)

bin_limits = np.arange(CLIP_LOWER, CLIP_UPPER+BIN_WIDTH, BIN_WIDTH)

n, bins, patches = plt.hist(residuals, bin_limits, facecolor='green', alpha=0.5)

print("         Bin             Count")
print("-----------------------------------")
for bin_num in range(0, len(bin_limits)-1):
    percentage = 100 *(n[bin_num] / count)
    print("{:5.1f}dB <--> {:5.1f}dB: {:4.0f} ({:.1f}%)".format(bins[bin_num], bins[bin_num+1], n[bin_num], percentage))
print("\nMean: {:.2f} SD: {:.2f} RMSE: {:.2f}".format(mean, std_dev, rmse))

plt.ylim([0,2500])

plt.axvline(x=mean, color='r')
plt.axvline(x=mean-std_dev, ls='dashed', color='r')
plt.axvline(x=mean+std_dev, ls='dashed', color='r')

plt.xlabel('Error (dB)'.format(len(values), mean, std_dev))
plt.ylabel('Count')
#plt.title('Predicted vs. Measured Field Strength')
bbox_props = dict(fc="white", ec="k", lw=1)
summary = "Count = {:d}\nMean = {:.2f}\n$\sigma$ = {:.2f}".format(count, mean, std_dev)
plt.text(-55,2000, summary, bbox=bbox_props)
plt.show()
