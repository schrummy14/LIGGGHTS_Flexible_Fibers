# LIGGGHTS Flexible Fibers

[Home](Home)

[Commands](commands)

## Bond Style Harmonic Command

### Syntax

```text
bond_type harmonic
```

### Examples

```text
bond_type harmonic
bond_coeff 1 80.0 1.2
```

### Description

The *harmonic* bond style uses the following force and moment equation

![bond equation](equations/bond_harmonic.jpg "Equation")

Where *r0* is the equilibrium bond distance. Note that the usual 1/2 factor is included in *K*.

The following coefficeints must be defined for each bond type via the [bond_coeff](bond_coeff) command as in the example above, or in the data file or restart files by the [read_data](not_done_yet) or [read_restart](not_done_yet) commands:

* K (energy/distance^2)
* r0 (distance)

## Restrictions

This bond style can only be used if LIGGGHTS(R)-PUBLIC was built with the MOLECULAR package (which it is by default). See the [How to Install](how_to_install) section for more info on packages.

## Related Commands

[bond_coeff](bond_coeff)
[delete_bonds](not_done_yet)

### Defaults

none
