# LIGGGHTS Flexible Fibers

[Home](Home)

[Commands](commands)

## Bond Style Command

### Syntax

```text
bond_style style args
```

* style = none for any style except *hybrid*
* *hybrid* args = list of one or more styles

### Examples

```text
bond_style gran
bond_style harmonic
bond_style hybrid harmonic fene
```

### Description

Set the formula(s) LIGGGHTS(R)-PUBLIC uses to compute bond interactions
between pairs of atoms. In LIGGGHTS(R)-PUBLIC, a bond differs from a pairwise
interaction, which are set via the [pair_style](pair_style) command. Bonds
are defined between specified pairs of atoms and remain in force for the
duration of the simulation (unless the bond breaks which is possible in some
bond potentials). The list of bonded atoms is read in by a
[read_data](read_data) or [read_restart](read_restart) command from a data
or restart file. By contrast, pair potentials are typically defined between
all pairs of atoms within a cutoff distance and the set of active interactions
changes over time.

Hybrid models where bonds are computed using different bond potentials can be
setup using the *hybrid* bond style.

The coefficients associated with a bond style can be specified in a data or
restart file or via the [bond_coeff](bond_coeff) command.

All bond potentials store their coefficient data in binary restart files which
means bond_style and [bond_coeff](bond_coeff) commands do not need to be
re-specified in an input script that restarts a simulation. See the
[read_restart](read_restart) command for details on how to do this. The one
exception is that bond_style hybrid only stores the list of sub-styles in the
restart file; bond coefficients need to be re-specified.

---

#### Warning

When both a bond and pair style is defined, the [special_bonds](not_done_yet)
command often needs to be used to turn off (or weight) the pairwise
interaction that would otherwise exist between 2 bonded atoms.

### Restrictions

Bond styles can only be set for atom styles that allow bonds to be defined.

Most bond styles are part of the MOLECULAR package. They are only enabled if
LIGGGHTS(R)-PUBLIC was built with that package. See the
[Making LIGGGHTS(R)-PUBLIC](how_to_install) section for more info on packages.
The doc pages for individual bond potentials tell if it is part of a package.

### Related Commands

[gran](bond_gran),
[harmonic](bond_harmonic),
[hybrid](bond_hybrid),
[bond_coeff](bond_coeff),
[delete_bonds](not_done_yet)

### Default

bond_style none

© Copyright 2016, DCS Computing GmbH, JKU Linz and Sandia Corporation
