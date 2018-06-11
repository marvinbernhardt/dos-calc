# dos-calc
Calculation of translational, rotational and vibrational density of states

## Usage
The program expects parameters on stdin. Use -v to get informed, which parameter is read next.
The parameters are as follows:
nsamples
nblocks
nblocksteps
nmoltypes
moltype1_nmols moltype1_natomspermol moltype1_atom1masses moltype1_rot_treat moltype1_abc_indicators
...

Add more lines if you have multiple moleculetypes. 

## Limitations
PBC recombination does not work for non orthorhombic boxes. However you can unwrap the trajectory before using dos-calc and then use the --no-pbc option.
Rotation limitations ...
