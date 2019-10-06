# DosCalc

Calculation of translational, rotational and vibrational density of states

## Dependencies

- CBLAS
- LAPACKE
- FFTW
- [libxdrfile](https://github.com/wesbarnett/libxdrfile) 
- [cJSON](https://github.com/DaveGamble/cJSON) 

## Installation

```bash
# needed only if libxdrfile and/or cJSON are not installed globally
CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}:/path/to/libxdrfile 
CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}:/path/to/cJSON
export CMAKE_PREFIX_PATH
INSTALL_PREFIX=/destination/for/dos-calc

mkdir build
pushd build
cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX ..
make
make install
popd
```

By default RPATH is used for the location of libxdrfile, but you can turn it of using `-DUSE_RPATH=OFF`.
If turned off, libxdrfile.so and libcjson.so.1 have to be in a standard directory or in `LD_LIBRARY_PATH`.

When cmake is run with `-DDEBUG=ON` a lot of intermediate results will be printed during runtime.
This will mostly be useful for small test systems.

There is also a scripts folder, but those scripts are not automatically installed anywhere.

## Usage

Check `dos-calc --help` for command line options.
In any case you need to provide a parameter file in JSON format (e.g. `params.json`) and a trajectory in trr format (e.g. `traj.trr`).
The output will be written to the JSON file `dos.json` or whatever filename specified with `-o`.

## Parameter JSON

A short example:
```
{
    "nsamples": 5,
    "nblocks": 1,
    "nblocksteps": 1000,
    "moltypes": [
        {
            "nmols": 1950,
            "atom_masses": [ 15.9994, 1.008, 1.008],
            "rot_treat": "f",
            "abc_indicators": [1, 2, 0, -1]
        },
        {
            "nmols": 50,
            "atom_masses": [ 15.9994, 1.008, 15.035 ],
            "rot_treat": "f",
            "abc_indicators": [1, 2, 0, -1]
        }
    ],
    "cross_spectra": [
        {
            "name": "water_trn water_vib",
            "type": "i",
            "dof_pairs":
            [
                [{"moltype": 0, "dof_type": "t", "dof_list": [0, 1, 2]}, {"moltype": 0, "dof_type": "v", "dof_list": [0, 1, 2, 3, 4, 5, 6, 7, 8]}]
            ]
        },
    ]
}
```

The order of moltypes and the atoms therein must reflect your trajectory!

The `type` of a `cross spectrum` can be either "i" for inside molecule.
This will result in a cross correlation of degrees of freedom within molecules and an average over molecules.
The other option is "e" for external, which will correlate all specified degrees of freedom from all matching molecules (can take a long time).

The variable `dof_type` can be one of:

- "t" for translational xyz-component of the molecule. Possible indices for dof_list: 0, 1, 2.
- "r" for rotational xyz-component of each atom in the molecule. Possible indices for dof_list: 0, 1, ..., 3*mol_natoms.
- "v" for vibrational xyz-component of each atom in the molecule. Possible indices for dof_list: 0, 1, ..., 3*mol_natoms.
- "o" for rotational abc-component of the molecule. Possible indices for dof_list: 0, 1, 2.

Here, mol_natoms indicates the number of atoms in one molecule of the chosen moleculetype.


## Samples and Blocks

The following visualisation shows you the effect of `nsamples = 2` and `nblocks = 3`.

```
---------------------trajectory----------------------
----------sample0--------- ----------sample1---------
-block0- -block1- -block2- -block0- -block1- -block2-
```
In the consequence the trajectory must have equal or more than `nsamples * nblocks * nblocksteps` frames with positions and velocities.

For each sample DosCalc will generate DoS files. Each sample can consist of multipe blocks that contribute to the sample's DoS (for example to reduce noise).

## Rotational treatment

Every molecule has three principal axes of rotation.
They are calculated in DosCalc but the algorithm depends on the moleules point group.
You set two variables for this:
- `moltypeX_rot_treat`, which is a char that determines the method used
- `moltypeX_abc_indicators`, which is four int numbers, that determine two pairs of atoms (zero indexed).
  If the number is -1 it stands for the center of mass of the molecule.
  These atom pairs define the helping vectors a, b and c, that are related to the true principal axes (from lowest to highest moment of inertia)
  The first two numbers set a, the second two set b'. c is the cross product of a and b', b is the cross product of c and a.

- for atoms and monoatomic ions `moltypeX_rot_treat` and `moltypeX_abc_indicators` both will be ignored, if `moltypeX_natomspermol == 1`.
- for linear molecules set `moltypeX_rot_treat = 'l'` and `moltypeX_abc_indicators = 0 0 0 0` (the latter will be ignored).
  The rotational DoS will be with regard to the axes x, y, and z.
  Tested only for diatomic molecules.
- for molecules where the axis can swap order by vibration, for example ammonia, set `moltypeX_rot_treat = 'a'` and `moltypeX_abc_indicators = 1 2 0 -1`.
  Ammonia has atoms N H H H. Therefore 1 2 defines vector a between two hydrogens. 0 -1 defines vector b' from the nitrogen atom to the COM and thereby along the symmetry axis.
  You can not use the actual principal axes of rotation of the molecule, because they swap order during vibration.
  The rotational DoS will be with regard to a, b and c.
  Note, that this will yield unusable results, if a, b and c are not close to the actual principal axes.
- for other molecules, for example water, set `moltypeX_rot_treat = 'f'` and `moltypeX_abc_indicators = 1 2 0 -1`.
  Water has atoms O H H. Therefore 1 2 defines vector a between the two hydrogens. 0 -1 defines vector b' along the symmetry axis.
  The rotational DoS will be with regard the actual principal axes of rotation.
  a, b and c will be used to check, that the actual axis derived from the moment of inertia tensor, always point in the same direction.
- if you do not want any velocity separation use `moltypeX_rot_treat = 'u'` and the unseparated DoS will be output to DoS_vib.

## Limitations

- PBC recombination does not work for non orthorhombic boxes. However you can make all molecules whole in the trajectory before using dos-calc and then use the --no-pbc option.
- The program uses single precision.
- Linear molecules are assumed to be diatomic.
- Not tested on periodic molecules.

## Scripts

### show-doses.py

Check `show-doses.py --help` for command line options.
