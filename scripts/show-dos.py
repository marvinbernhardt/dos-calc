#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
from pathlib import Path
import sys

if sys.version_info < (3, 6):
    print("This script requires at least Python version 3.6")
    sys.exit(1)

parser = argparse.ArgumentParser(description='Plot DoS files.')
parser.add_argument('-c', '--components', help='show components DoS',
                    dest='show_components', action='store_true')
parser.add_argument('-x', '--cross', help='show cross DoS',
                    dest='show_cross', action='store_true')
parser.add_argument('-s', '--samples', help='show samples individually',
                    dest='show_samples', action='store_true')

args = parser.parse_args()

freq_file = Path(f"frequencies.txt")
if freq_file.is_file():
    frequencies = np.array(pd.read_csv(f"frequencies.txt", sep=' ', header=None)).flat
    print("frequencies found")
else:
    frequencies = False
    print("no frequencies found")

doses = [
    {'name': 'dos_trn', 'tags': ['standard']},
    {'name': 'dos_rot', 'tags': ['standard']},
    {'name': 'dos_vib', 'tags': ['standard']},
    {'name': 'dos_trn_x', 'tags': ['comp']},
    {'name': 'dos_trn_y', 'tags': ['comp']},
    {'name': 'dos_trn_z', 'tags': ['comp']},
    {'name': 'dos_rot_a', 'tags': ['comp']},
    {'name': 'dos_rot_b', 'tags': ['comp']},
    {'name': 'dos_rot_c', 'tags': ['comp']},
    {'name': 'dos_vib_x', 'tags': ['comp']},
    {'name': 'dos_vib_y', 'tags': ['comp']},
    {'name': 'dos_vib_z', 'tags': ['comp']},
    {'name': 'dos_x_trn_rot', 'tags': ['cross']},
    {'name': 'dos_x_trn_vib', 'tags': ['cross']},
    {'name': 'dos_x_rot_vib', 'tags': ['cross']},
]

samples = [int(dos_trn_file.split('_')[0].split('e')[1]) for dos_trn_file in glob.glob('sample*_dos_trn.txt')]

if len(samples) == 0:
    raise Exception("ERROR: No dos found to show!")
elif len(samples) == 1 or args.show_samples:
    for sample in samples:

        # filter which doses to show
        doses_to_show = []
        for dos in doses:
            if 'standard' in dos['tags']:
                doses_to_show.append(dos)
            if 'comp' in dos['tags'] and args.show_components:
                doses_to_show.append(dos)
            if 'cross' in dos['tags'] and args.show_cross:
                doses_to_show.append(dos)

        # load data from files
        for dos in doses_to_show:
            dos['data'] = np.array(pd.read_csv(f"sample{sample}_{dos['name']}.txt", sep=' ', header=None))

        # show dos integrals
        for moltype in range(len(doses_to_show[0]['data'])):
            if frequencies:
                print(f"sample {sample}, moltype {moltype}")
            else:
                factor = 1 / (np.trapz(doses_to_show[0]['data'][moltype]) / 3)
                print(f"sample {sample}, moltype {moltype}, normalized to dos_trn/3")

            for dos in doses_to_show:
                if frequencies:
                    print(f"integral {dos['name']}:", np.trapz(dos['data'][moltype], frequencies))
                else:
                    print(f"integral {dos['name']}:", np.trapz(dos['data'][moltype]) * factor)

        # show dos plots
        fig, ax = plt.subplots()
        ax.set_title(f"sample {sample}")
        for moltype in range(len(doses_to_show[0]['data'])):
            for dos in doses_to_show:
                if frequencies:
                    ax.plot(frequencies, dos['data'][moltype], label=f"{moltype} {dos['name']}")
                else:
                    ax.plot(dos['data'][moltype], label=f"{moltype} {dos['name']}")

        if frequencies:
            plt.xlabel("freqency in 1/ps")
        fig.legend()

    plt.show()

else:
    # filter which doses to show
    doses_to_show = []
    for dos in doses:
        if 'standard' in dos['tags']:
            doses_to_show.append(dos)
        if 'comp' in dos['tags'] and args.show_components:
            doses_to_show.append(dos)
        if 'cross' in dos['tags'] and args.show_cross:
            doses_to_show.append(dos)

    # load data from files
    for dos in doses_to_show:

        dos['data'] = []
        for sample in samples:
            dos['data'].append(np.array(pd.read_csv(f"sample{sample}_{dos['name']}.txt", sep=' ', header=None)))

        dos['mean'] = np.array(dos['data']).mean(axis=0)
        dos['std'] = np.array(dos['data']).std(axis=0)

    # show dos integrals
    for moltype in range(len(doses_to_show[0]['data'][0])):
        if frequencies:
            print(f"moltype {moltype}")
        else:
            factor = 1 / (np.trapz(doses_to_show[0]['mean'][moltype]) / 3)
            print(f"moltype {moltype}")
            print(f"sample {sample}, integrals normalized to dos_trn/3")

        for dos in doses_to_show:
            if frequencies:
                print(f"integral {dos['name']}:", np.trapz(dos['mean'][moltype], frequencies))
            else:
                print(f"integral {dos['name']}:", np.trapz(dos['mean'][moltype]) * factor)

    # show dos plots
    fig, ax = plt.subplots()
    for moltype in range(len(doses_to_show[0]['data'][0])):
        for dos in doses_to_show:
            if frequencies:
                line, = ax.plot(frequencies, dos['mean'][moltype], label=f"{moltype} {dos['name']}")
                ax.fill_between(frequencies,
                                dos['mean'][moltype]-dos['std'][moltype],
                                dos['mean'][moltype]+dos['std'][moltype],
                                facecolor=line.get_color(), alpha=0.3)
            else:
                arbitrary_freqs = np.arange(0, len(dos['mean'][moltype]))
                line, = ax.plot(dos['mean'][moltype], label=f"{moltype} {dos['name']}")
                ax.fill_between(arbitrary_freqs,
                                dos['mean'][moltype]-dos['std'][moltype],
                                dos['mean'][moltype]+dos['std'][moltype],
                                facecolor=line.get_color(), alpha=0.3)

    if frequencies:
        plt.xlabel("freqency in 1/ps")
        plt.ylabel("DoS in ps")
    else:
        plt.xlabel("freqency in arbitrary units")
        plt.ylabel("DoS in arbitrary units")
    fig.legend()
    plt.show()
