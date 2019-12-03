#!/bin/bash
rm restart_*
rm post/*

echo "Bash Version = ${BASH_VERSION}"
for i in {1..10}
do
	mpirun -n $i liggghts_bonds_3_8 -in in.liggghts
	mv post post$i
done
