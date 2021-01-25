# LIGGGHTS Flexible Fibers

[Home](Home)

[Commands](commands)

## Fix Void Ratio Command

### Syntax

```text
fix ID group-ID voidratio region regArg ntry ntryArg nevery neveryArg
```

* ID, group-ID are documented in the [fix](fix) command
* voidratio = style name of this fix command
* region regArg = region name for the void ratio to be calculated
* ntry ntryArg = the number of random points used to calculate the void ratio
* nevery neveryArg = the number of steps between void ratio calculations (set
    to 0 for the calculation to be done at the end of a run)

### Examples

```text
fix vR all voidratio region reg ntry 1000000 nevery 5000
```

### Description

Calculates the void ratio inside a specific region durring the simulation. The
void ratio can be accessed via

```text
variable voidRatio equal f_vR[1]
```

The void ratio is calculated locally on each individual thread and then summed
together to give the total void ratio. This is done by randomly selecting a
point within the processor's domain and checking if that point is inside an
atom.

### Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to
[binary restart files](restart). None of the [fix_modify](fix_modify)
options are relevant to this fix. No global or per-atom quantities are stored
by this fix for access by various [output commands](dump). No
parameter of this fix can be used with the start/stop keywords of the
[run](run) command.

### Restrictions

None

#### Defaults

None
