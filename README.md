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

By default RPATH is used for the location of libxdrfile.so and libcjson.so.1, but it can be turned off using `-DUSE_RPATH=OFF`.
If turned off, libxdrfile.so and libcjson.so.1 have to be in a standard directory or in `LD_LIBRARY_PATH` at runtime.

When cmake is run with `-DDEBUG=ON` a lot of intermediate results are printed during runtime.
This is mostly be useful for very small test systems.

There is also a scripts folder, but those scripts are not automatically installed anywhere.

## Usage

Check `dos-calc --help` for command line options.
In any case one needs to provide a parameter file in JSON format (e.g. `params.json`) and a trajectory in trr format (e.g. `traj.trr`).
The output is be written to the JSON file `dos.json` or whatever filename specified with `-o`.

## Input

A example params.json of a mixture of three-point model water with united atom methanol.
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
        }
    ]
}
```

### Samples and Blocks

The following visualisation shows the effect of `nsamples = 2` and `nblocks = 3`.

```
---------------------trajectory----------------------
----------sample0--------- ----------sample1---------
-block0- -block1- -block2- -block0- -block1- -block2-
```
In the consequence the trajectory must have equal or more than `nsamples * nblocks * nblocksteps` frames with positions and velocities.

For each sample DosCalc is generating power spectra in the output file. Each sample can consist of multipe blocks that contribute to the sample's DoS (for example to reduce noise).

### Moltypes (Topology)

The order of moltypes and the atoms must reflect the structure of the trajectory!

- `nmols` is the number of molecules of that type.
- `atom_masses` indicates the atoms of each molecule of that type.

One can define less atoms in total, than are present in the trajectory (produces a warning) but not more (produces an error).

### Rotational treatment

For each moltype there is also:

- `rot_treat`, which is a char that determines the method used for the rotational treatment.
- `abc_indicators`, which is a list of four integers, which indicate two pairs of atoms (zero indexed).
  If the number is -1 it stands for the center of mass of the molecule.
  These atom pairs define the helping vectors a, b and c, that are related to the true principal axes (from lowest to highest moment of inertia)
  The first two numbers define a, the second two define b'. c is the cross product of a and b', b is the cross product of c and a.

Depending on the point group of the molecule the following needs to be defined:

- For atoms and monoatomic ions (if `atom_masses` has only one element) `rot_treat` and `abc_indicators` both are ignored.
- For linear molecules set `"rot_treat": "l"`. `abc_indicators` is ignored.
  The rotational DoS is calculated with regard to the axes x, y, and z (but still be called `rot_a`, ... in the output file.
  Tested only for diatomic molecules.
- For molecules where the axis can swap order by vibration, for example ammonia, set `"rot_treat": "a"` and `"abc_indicators": [1, 2, 0 -1]`.
  Ammonia has atoms N H H H. Therefore 1 2 defines vector a between two hydrogens. 0 -1 defines vector b' from the nitrogen atom to the COM and thereby along the symmetry axis.
  One can not use the actual principal axes of rotation of the molecule, because they swap order during vibration.
  The rotational DoS is calculated be with regard to a, b and c.
  Note, that this does yield unusable results, if a, b and c are not close to the actual principal axes.
- For other molecules, for example water, set `"rot_treat": "f"` and `"abc_indicators": [1, 2, 0 -1]`.
  Water has atoms O H H. Therefore 1 2 defines vector a between the two hydrogens. 0 -1 defines vector b' along the symmetry axis.
  The rotational DoS is calculated with regard to the actual principal axes of rotation.
  The vectors a, b and c are used to ensure, that the actual axis derived from the moment of inertia tensor, always point in the same direction.
  This will likely be the correct choice for most molecules, especially larger ones.
- If one does not want any velocity separation `"rot_treat": "u"` can be used and the unseparated DoS is written to `vib_{x,y,z}`.

### Cross spectrum

The `name` of a `cross spectrum` can be up to 79 characters long and can help identify the specified cross correlation in the output file.
The `type` of a `cross spectrum` can be "i" for inside or "e" for external.
Option "i" does result in a cross correlation of degrees of freedom within molecules and an average over molecules.
Option "e" does correlate all specified degrees of freedom from all matching molecules (can take a long time).

The variable `dof_type` can be one of:

- "t" for translational xyz-component of the molecule. Possible indices for dof_list: 0, 1, 2.
- "r" for rotational xyz-component of each atom in the molecule. Possible indices for dof_list: 0, 1, ..., 3*mol_natoms.
- "v" for vibrational xyz-component of each atom in the molecule. Possible indices for dof_list: 0, 1, ..., 3*mol_natoms.
- "o" for rotational abc-component of the molecule. Possible indices for dof_list: 0, 1, 2.

Here, mol_natoms indicates the number of atoms in one molecule of the chosen moleculetype.

## Output

A example dos.json corresponding to the example above with long lists of numbers shortened to [...].
```
{
    "frequencies":	[...],
    "moltypes":	[{
        "dos":	[{
            "name":	"trn_x",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"trn_y",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"trn_z",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"rot_x",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"rot_y",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"rot_z",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"vib_x",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"vib_y",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"vib_z",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"rot_a",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"rot_b",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"rot_c",
            "data":	[[...], [...], [...], [...], [...]]
        }],
        "moments_of_inertia":	[[0.0071233315393328667, 0.013312174007296562, 0.020435504615306854], [0.0071259848773479462, 0.013308963738381863, 0.02043495886027813], [0.0071264943107962608, 0.013315362855792046, 0.020441846922039986], [0.0071183349937200546, 0.013314731419086456, 0.020433064550161362], [0.007123015820980072, 0.013312199153006077, 0.020435240119695663]]
    }, {
        "dos":	[{
            "name":	"trn_x",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"trn_y",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"trn_z",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"rot_x",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"rot_y",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"rot_z",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"vib_x",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"vib_y",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"vib_z",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"rot_a",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"rot_b",
            "data":	[[...], [...], [...], [...], [...]]
        }, {
            "name":	"rot_c",
            "data":	[[...], [...], [...], [...], [...]]
        }],
        "moments_of_inertia":	[[0.00785850640386343, 0.16801847517490387, 0.17587696015834808], [0.0078581739217042923, 0.16808962821960449, 0.17594780027866364], [0.0078575713559985161, 0.16799667477607727, 0.17585422098636627], [0.00785425677895546, 0.16799832880496979, 0.17585258185863495], [0.007849549874663353, 0.16799871623516083, 0.17584823071956635]]
    }],
    "cross_spectra":	[{
        "name":	"water_trn water_vib",
        "data":	[[...], [...], [...], [...], [...]]
    }]
}
```

The first dimension of each `data` and `moments_of_inertia` list is determined by `nsamples`.
The second dimension of each `data` list is the number of frequencies, which is `floor(nblocksteps / 2) + 1`.

A verbal (and therefore not exact) description of the DoS components:

- `trn_x` is the power spectrum of the x component of the velocities of the molecules center of mass translational motion.
- `rot_x` is the power spectrum of the x component of the velocities of the atoms due to molecular rotational motion.
- `vib_x` is the power spectrum of the x component of the velocities of the atoms due to molecular vibrational motion.
- `rot_a` is the power spectrum of the angular velocity around principal axis 'a' of each molecule.

## Limitations

- PBC recombination does not work for non-orthorhombic boxes.
  However one can make all molecules whole in the trajectory (`gmx trjconv -pbc mol`) before using `dos-calc` and then use the --no-pbc option.
  Whole molecules are alowed to jump between frames so a full unwrapping of the trajectory is not necessary.
- The program uses single precision.
- Linear molecules are assumed to be diatomic.
- Not tested and likely to fail on periodic molecules.

## Scripts

### show-doses.py

Check `show-doses.py --help` for command line options.
