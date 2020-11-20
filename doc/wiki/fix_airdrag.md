# LIGGGHTS Flexible Fibers

[Home](Home)

[Commands](commands)

## Fix Air Drag Command

### Syntax

```text
fix ID group-ID airdrag fluid_viscocity fluid_density wx wy wz region
```

* ID, group-ID are documented in the [fix](fix) command
* airdrag = style name of this fix command
* fluid_viscocity = viscocity of fluid (Pa s)
* fluid_density = density of fluid (kg/m3)

*Optional* *Inputs*

* wx, wy, wz = fluid velocity in the x, y, and z directions (m/s)
* region = region name for the fix to interact with

### Examples

```text
fix 1 all airdrag 1.8e-5 1.23
fix 1 all airdrag 1.8e-5 1.23 0.0 0.0 1.0
fix 1 all airdrag 1.8e-5 1.23 0.0 0.0 1.0 vertChamber
```

### Description

Add a fluid to the simulation. The fluid interacts with the particles by the following equations outlined by "Classical Mechanics" by John R. Taylor for forces and "Viscous torque on a sphere under arbitrary rotation" by U. Lei er al. for torque calculations.

![air drag equation](equations/fix_airdrag.jpg "Equation")

---

#### Warning

The units for the coefficients are currently static and the user must make sure that they are following what is laid out in this document. No options for variable input are allowed at this time (v_wx). Currently to change flow properties, one must unfix and then apply a new fix.

This fix is NOT a subsitute for a true CFD simulation.

---

### Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart files](not_done_yet). None of the [fix_modify](not_done_yet) options are relevant to this fix. No global or per-atom quantities are stored by this fix for access by various [output commands](not_done_yet). No parameter of this fix can be used with the start/stop keywords of the [run](not_done_yet) command.

The forces due to this fix are imposed during an energy minimization, invoked by the minimize command. This fix should only be used with damped dynamics minimizers that allow for non-conservative forces. See the min_style command for details.

### Restrictions

None

### Related Commands

[fix viscous](not_done_yet)

#### Defaults

wx = wy = wz = 0.0,
region = simulation domain