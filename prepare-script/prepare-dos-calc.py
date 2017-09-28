#!/usr/bin/env python3

import gromacstools as gt

# modify this stuff to your needs!
atomtypes = [{'name': 'C',  'mass': 12.011},
             {'name': 'H',  'mass':  1.008},
             {'name': 'CL', 'mass': 35.453},
             {'name': 'CH', 'mass': 13.019}]

simpletop = [{'name': "CHC", 'atomtypes': [3, 2, 2, 2], 'nmols': 512,
              'sigma': 3, 'abc_indicators': [2, 3, 1, -1]}]

nsamples = 1
nblocks = 5
nblocksteps = 5000
outputfile = "params.txt"
# stop modifying here

top = gt.top.Topology()
top.load_simple_top(simpletop, atomtypes)
top.save_2pt_parameters_file(nsamples, nblocks, nblocksteps,
                             parameters_file=outputfile)
