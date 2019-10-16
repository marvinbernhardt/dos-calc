#!/usr/bin/env python3

import argparse
import json
import matplotlib.pyplot as plt
import numpy as np
import sys

if sys.version_info < (3, 6):
    print("This script requires at least Python version 3.6")
    sys.exit(1)

    
def create_nocomp_spectra(system):
    # definitions
    dos_nocomp_dict = {
        'trn': ['trn_x', 'trn_y', 'trn_z'],
        'rot': ['rot_x', 'rot_y', 'rot_z'],
        'vib': ['vib_x', 'vib_y', 'vib_z'],
        'roto': ['roto_a', 'roto_b', 'roto_c']
    }

    # sum no-component spectra from component spectra
    for h, moltype in enumerate(system['moltypes']):
        # loop dos_types to build
        for dos_type_nocomp, dos_types_comp in dos_nocomp_dict.items():
            dos_nocomp_data = None
            # loop dos_types to sum
            for dos_type_comp in dos_types_comp:
                dos_comp_data = np.array(moltype['spectra'][dos_type_comp])
                if dos_nocomp_data is None:
                    dos_nocomp_data = dos_comp_data.copy()
                else:
                    dos_nocomp_data += dos_comp_data
            moltype['spectra'][dos_type_nocomp] = dos_nocomp_data
        

def _plot_spectra(system, show_components, show_cross, show_samples, show_roto, temperature, sample=None):

    K_GRO = 0.00831445986144858
    # frequencies
    freq = system['frequencies']

    # sizes
    nmoltypes = len(system['moltypes'])

    linestyles = ['--', '-.', ':',
                  (0, (3, 2, 1, 2, 1, 2))
                 ]

    fig, ax = plt.subplots(figsize=(5, 3))
    for h, moltype in enumerate(system['moltypes']):
        # filter spectra for comp/nocomp and roto
        filtered_spectra = {dos_name:dos for (dos_name, dos) in moltype['spectra'].items()
                            if not (not show_components and '_' in dos_name)
                            if not (show_components and not '_' in dos_name)
                            if not (not show_roto and dos_name.startswith('roto'))}
        nspectra = len(filtered_spectra)
        # loop dos_types to plot
        for d, (dos_name, dos) in enumerate(filtered_spectra.items()):
            if show_samples:
                dos = np.array(dos)[([sample], slice(None))]
            else:
                dos = np.array(dos)
            # show integral
            if temperature is None:
                prefactor = 1
                print(f"integral {dos_name}: {np.trapz(np.mean(dos, axis=0), freq):.4E} kJ/mol")
            else:
                prefactor = 1 / K_GRO / temperature
                print(f"integral {dos_name}: {np.trapz(np.mean(dos, axis=0), freq) * 2 / K_GRO / temperature:.4E}")
            # plot
            if show_components:
                color = plt.get_cmap('brg')(h/nmoltypes+d/nmoltypes/nspectra)
                linestyle=linestyles[d%3]
            else:
                color = plt.get_cmap('Dark2')(h/8)
                linestyle=linestyles[d%len(linestyles)]
            line, = ax.plot(freq, prefactor * np.mean(dos, axis=0),
                            linestyle=linestyle,
                            color=color,
                            label=f"{h} {dos_name}")
            if not show_samples:
                ax.fill_between(freq,
                                prefactor * np.min(dos, axis=0),
                                prefactor * np.max(dos, axis=0),
                                facecolor=line.get_color(), alpha=0.3)

    if show_cross:
        for d, (cs_name, cs) in enumerate(system['cross_spectra'].items()):
            if show_samples:
                cs = np.array(cs)[([sample], slice(None))]
            else:
                cs = np.array(cs)
            # plot
            if temperature is None:
                prefactor = 1
                print(f"integral {cs_name}: {np.trapz(np.mean(cs, axis=0), freq):.4E} kJ/mol")
            else:
                prefactor = 1 / K_GRO / temperature
                print(f"integral {cs_name}: {np.trapz(np.mean(cs, axis=0), freq) * 2 / K_GRO / temperature:.4E}")
            line, = ax.plot(freq, prefactor * np.mean(cs, axis=0),
                            color=plt.get_cmap('Set2')(d/8),
                            label=cs_name)
            if not show_samples:
                ax.fill_between(freq,
                                prefactor * np.min(cs, axis=0),
                                prefactor * np.max(cs, axis=0),
                                facecolor=line.get_color(), alpha=0.3)
            
    ax.set_xlabel("freqency in 1/ps")
    if temperature is None:
        ax.set_ylabel("DoS in kJ/mol ps")
    else:
        ax.set_ylabel("DoS in ps")
    ax.set_xlim(0)
    ax.set_ylim(0)
    ax.legend()


def plot_spectra(system, show_components, show_cross, show_samples, show_roto, temperature):
    nsamples = len(system['moltypes'][0]['spectra']['trn_x'])

    if nsamples == 1 or show_samples:
        for sample in range(nsamples):
            _plot_spectra(system, show_components, show_cross, show_samples, show_roto, temperature, sample)
        plt.show()
    else:
        _plot_spectra(system, show_components, show_cross, show_samples, show_roto, temperature)
        plt.show()

if __name__ == "__main__":

    # arguemnt parser
    parser = argparse.ArgumentParser(description='Plot DoS files.')
    parser.add_argument('filename', nargs='?', help='file to plot DoS from',
                        type=argparse.FileType('r'), default='dos.json')
    parser.add_argument('-c', '--components', help='show components DoS',
                        dest='show_components', action='store_true')
    parser.add_argument('-x', '--cross', help='show cross DoS',
                        dest='show_cross', action='store_true')
    parser.add_argument('-s', '--samples', help='show samples individually',
                        dest='show_samples', action='store_true')
    parser.add_argument('-o', '--roto', help='show DoS_rotomega',
                        dest='show_roto', action='store_true')
    parser.add_argument('-t', '--temperature', help='show integral normalized by 1/(kT)',
                        dest='temperature', type=float, default=None)

    args = parser.parse_args()
    system = json.load(args.filename)
    create_nocomp_spectra(system)
    plot_spectra(system, args.show_components, args.show_cross, args.show_samples, args.show_roto, args.temperature)
