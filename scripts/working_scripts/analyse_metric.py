import os
from os import path
import sys
import json
import requests
import argparse
from collections import defaultdict
from typing import Tuple
import numpy as np
from matplotlib import pyplot as plt

#  CONSTANTS  ##############################################################################

from constants import *

#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('metric_tsv', help='TSV file with columns label, metric_value')
args = parser.parse_args()

metric_tsv = args.metric_tsv

#  MAIN  ##############################################################################

labels = []
values = []
with open(metric_tsv) as r:
    for line in r:
        domain, label, value = line.split()
        labels.append(label)
        values.append(float(value))

labels = np.array(labels)
values = np.array(values)

bins = np.arange(0, 30, 0.5)

plt.hist(values, bins=bins)
plt.title('All SSEs')
# plt.show()
plt.savefig(path.join('metric_histograms', 'all.png'))
plt.clf()

print('label', '<10', '<20', sep='\t')
for label in sorted(set(labels)):
    vals = values[labels == label]
    plt.hist(vals, bins=bins)
    plt.title(label)
    # plt.show()
    plt.savefig(path.join('metric_histograms', label + '.png'))
    plt.clf()
    print(label, (vals <= 10).sum() / vals.shape[0], (vals <= 20).sum() / vals.shape[0], sep='\t')

for threshold in [5, 10, 12, 15, 20, 25, 30]:
    ratio = (values <= threshold).sum() / values.shape[0]
    print(f'<={threshold}:  {ratio}')
