# dos-calc

Calculation of translational, rotational and vibrational density of states

## Usage

Check `dos-calc --help` for command line options.
In any case you need to provide a trajectory with `-f`.

The program expects parameters on stdin. Use -v to get informed, which parameter is read next.
The parameters are as follows:
```
nsamples
nblocks
nblocksteps
nmoltypes
moltype1_nmols moltype1_natomspermol moltype1_atom1masses moltype1_rot_treat moltype1_abc_indicators
moltype2...
```

Add more lines if you have multiple moleculetypes.

I suggest to write the parameters in a file and pipe them into the program.

## Samples and Blocks

The following visualisation shows you the effect of `nsamples = 2` and `nblocks = 3`.

```
---------------------trajectory----------------------
----------sample0--------- ----------sample1---------
-block0- -block1- -block2- -block0- -block1- -block2-
```
For each sample dos-calc will generate DoS files. Each sample can consist of multipe blocks that contribute to the sample's DoS (for example to reduce noise).

## Rotational treatment

Every molecule has three principal axes of rotation.
They are calculated in dos-calc, but the algorithm depends on the moleules point group.
You set two variables for this:
* `moltypeX_rot_treat`, which is a char that determines the method used
* `moltypeX_abc_indicators`, which is four int numbers, that determine two pairs of atoms (zero indexed).
  If the number is -1 it stands for the center of mass of the molecule.
  These atom pairs define the helping vectors a, b and c, that are related to the true principal axes (from lowest to highest moment of inertia)
  The first two numbers set a, the second two set b'. c is the cross product of a and b', b is the cross product of c and a.

* for atoms and ions set `moltypeX_rot_treat = 'f'` and `moltypeX_abc_indicators = 0 0 0 0` (both will be ignored, if `moltypeX_natomspermol == 1`).
* for linear molecules set `moltypeX_rot_treat = 'l'` and `moltypeX_abc_indicators = 0 0 0 0` (the latter will be ignored).
  The rotational DoS will be with regard to the axes x, y, and z.
  Tested only for diatomic molecules.
* for molecules where the axis can swap order by vibration, for example ammonia, set `moltypeX_rot_treat = 'a'` and `moltypeX_abc_indicators = 1 2 0 -1`.
  Ammonia has atoms N H H H. Therefore 1 2 defines vector a between two hydrogens. 0 -1 defines vector b' from the nitrogen atom to the COM and thereby along the symmetry axis.
  The rotational DoS will be with regard to a, b and c.
  You can not use the actual principal axes of rotation of the molecule, because they swap order during vibration.
* for other molecules, for example water, set `moltypeX_rot_treat = 'f'` and `moltypeX_abc_indicators = 1 2 0 -1`.
  Water has atoms O H H. Therefore 1 2 defines vector a between the two hydrogens. 0 -1 defines vector b' along the symmetry axis.
  The rotational DoS will be with regard the actual principal axes of rotation.
  a, b and c will be used to check, that the actual axis derived from the moment of inertia tensor, always point in the same direction.

## Limitations

PBC recombination does not work for non orthorhombic boxes. However you can unwrap the trajectory before using dos-calc and then use the --no-pbc option.
The program uses single precision.
