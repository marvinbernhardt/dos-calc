#!/usr/bin/env python3.6

# modify these lines to your needs!
atomtypes = [{'name': 'O',  'mass': 15.9994},
             {'name': 'H',  'mass':  1.008}]

moltypes = [{'name': "SOL", 'atoms': [atomtypes[j] for j in [0, 1, 1]], 'nmols': 500,
             'sigma': 2, 'abc_indicators': [1, 2, 0, -1], 'rot_treat': 'f'},
            {'name': "H2O2", 'atoms': [atomtypes[j] for j in [1, 0, 0, 1]], 'nmols': 2,
             'sigma': 2, 'abc_indicators': [0, 3, 1, -1], 'rot_treat': 'f'}]

nsamples = 2
nblocks = 2
nblocksteps = 1000
outputfile = "params.txt"
# stop modifying here

with open(outputfile, 'w') as f:

    f.write(str(nsamples) + "\n")
    f.write(str(nblocks) + "\n")
    f.write(str(nblocksteps) + "\n")
    f.write(str(len(moltypes)) + "\n")
    for moltype in moltypes:
        moltype_atommasses = (atom['mass'] for atom in moltype['atoms'])
        moltype_string = f"{moltype['nmols']} "\
                         f"{len(moltype['atoms'])} "\
                         f"{' '.join(map(str, moltype_atommasses))} "\
                         f"{moltype['rot_treat']} "\
                         f"{' '.join(map(str, moltype['abc_indicators']))}"

        f.write(moltype_string + '\n')
