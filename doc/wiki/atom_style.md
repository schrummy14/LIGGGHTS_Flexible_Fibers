# LIGGGHTS Flexible Fibers

[Home](Home)

[Commands](commands)

## Atom Style Command

### Syntax

```text
atom_style style args
```

* style = *bond* or *charge* or *hybrid* or *sphere* or *granular* or *bond/gran* or *superquadric* or *sph*

```text
args = none for any style except hybrid and bond/gran
hybrid args = list of one or more sub-styles, each with their args
bond/gran args = n_bondtypes k and bonds_per_atom m and/or disableNormalContact n
    n_bondtypes k = number of bond types to be used in the simulation (k > 0)
    bonds_per_atom m = max number of bonds each atom can have in the simulation (m > 0)
    disableNormalContact n = disables the contact equations for a bonded pair of atoms if set to 1. (0 or 1)
```

### Examples

```text
* atom_style sphere

* atom_style superquadric

* atom_style hybrid sphere bond/gran n_bondtypes 1 bonds_per_atom 8 &
                           disableNormalContact 1
```

### Description

Define what style of atoms to use in a simulation. This determines what attributes are associated with the atoms. This command must be used before a simulation is setup via a [read_data](not_done_yet), [read_restart](not_done_yet), or [create_box](not_done_yet) command.

Once a style is assigned, it cannot be changed, so use a style general enough to encompass all attributes. E.g. with style bond, angular terms cannot be used or added later to the model. It is OK to use a style more general than needed, though it may be slightly inefficient.

The choice of style affects what quantities are stored by each atom, what quantities are communicated between processors to enable forces to be computed, and what quantities are listed in the data file read by the read_data command.

These are the additional attributes of each style and the typical kinds of physical systems they are used to model. All styles store coordinates, velocities, atom IDs and types. See the [read_data](not_done_yet), [create_atoms](not_done_yet), and [set](not_done_yet) commands for info on how to set these various quantities.

```text
bond          | bonds                            | bead-spring polymers
---------------------------------------------------------------------------
bond/gran     | number of bonds and bond info    | granular bond models
---------------------------------------------------------------------------
charge        | charge                           | atomic system with
              |                                  | charges
---------------------------------------------------------------------------
sph           | q(pressure), density             | SPH particles
---------------------------------------------------------------------------
sphere        | diameter, mass, angular velocity | granular models
---------------------------------------------------------------------------
granular      | diameter, mass, angular velocity | granular models
---------------------------------------------------------------------------
superquadric  | semi-axes, blockiness parameters,| granular models
              | mass, angular velocity,          |
              | quaternion                       |
---------------------------------------------------------------------------
```

#### Warning

It is possible to add some attributes, such as a molecule ID, to atom styles that do not have them via the [fix property/atom](not_done_yet) command. This command also allows new custom attributes consisting of extra integer or floating-point values to be added to atoms. See the [fix property/atom](not_done_yet) doc page for examples of cases where this is useful and details on how to initialize, access, and output the custom values.

---

All of the styles assign mass to particles on a per-type basis, using the [mass](not_done_yet) command, except for *sphere* or *granular* styles. They assign mass to individual particles on a per-particle basis.

For the *sphere* style, the particles are spheres and each stores a per-particle diameter and mass. If the diameter > 0.0, the particle is a finite-size sphere. If the diameter = 0.0, it is a point particle. This is typically used for granular models. Instead of sphere, keyword granular can be used.

For the *bond/gran* style, the number of granular bonds per atom is stored, and the information associated to it: the type of each bond, the ID of the bonded partner atom and the so-called bond history. The bond history is similar to the contact history for granular interaction, it stores the internal state of the bond. What exactly is stored in this internal state is defined by the granular [bond style](bond_style) used. There are 2 parameters: The number of bond types, and the maximum number of bonds that each atom can have. For each bond type, the parameters have to be specified via the [bond_coeff](bond_coeff) command. Note that [bond/gran](bond_gran) is an experimental code which is may not be available in your release of LIGGGHTS.

Typically, simulations require only a single (non-hybrid) atom style. If some atoms in the simulation do not have all the properties defined by a particular style, use the simplest style that defines all the needed properties by any atom. For example, if some atoms in a simulation are charged, but others are not, use the charge style. If some atoms have bonds, but others do not, use the bond style.

The only scenario where the hybrid style is needed is if there is no single style which defines all needed properties of all atoms. For example, if you want dipolar particles which will rotate due to torque, you would need to use “atom_style hybrid sphere dipole”. When a hybrid style is used, atoms store and communicate the union of all quantities implied by the individual styles.

LIGGGHTS(R)-PUBLIC can be extended with new atom styles as well as new body styles; see [this section](not_done_yet).

### Restrictions

This command cannot be used after the simulation box is defined by a [read_data](not_done_yet) or [create_box](not_done_yet) command.

The superquadric style is not yet available in the PUBLIC version The convexhull style is not yet available in the PUBLIC version

The bond, molecular styles are part of the MOLECULAR package. The line and tri styles are part of the ASPHERE package. They are only enabled if LIGGGHTS(R)-PUBLIC was built with that package. See the [Making LIGGGHTS(R)-PUBLIC](how_to_install) section for more info.

### Related Commands

[read_data](not_done_yet),
[pair_style](not_done_yet)

#### Default: none

##### © Copyright 2016, DCS Computing GmbH, JKU Linz and Sandia Corporation
