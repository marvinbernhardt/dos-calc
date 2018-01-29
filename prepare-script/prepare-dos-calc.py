#!/usr/bin/env python3

import gromacstools as gt

# modify this stuff to your needs!
atomtypes = [{'name': 'O',  'mass': 15.9994},
             {'name': 'H',  'mass':  1.008}]

simpletop = [{'name': "SOL", 'atoms': [atomtypes[j] for j in [0, 1, 1]], 'nmols': 5,
              'sigma': 2, 'abc_indicators': [1, 2, 0, -1], 'rot_treat': 'f'},
             {'name': "H2O2", 'atoms': [atomtypes[j] for j in [1, 0, 0, 1]], 'nmols': 2,
              'sigma': 2, 'abc_indicators': [0, 3, 1, -1], 'rot_treat': 'f'}]

nsamples = 2
nblocks = 2
nblocksteps = 1000
outputfile = "params.txt"
# stop modifying here

top = gt.top.Topology()
top.load_simple_top(simpletop)
top.save_2pt_parameters_file(nsamples, nblocks, nblocksteps,
                             parameters_file=outputfile)
