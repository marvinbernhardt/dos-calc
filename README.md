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
