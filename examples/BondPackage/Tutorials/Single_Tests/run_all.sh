
set -e

# Run all simulations

liggghts_prog=$1

num_cores=$2

cd cantilever_beam
echo 'Running cantiliever beam test'
rm -rf post/*
$liggghts_prog -in in.liggghts
cd ..

cd 3_point_bending
echo 'Running 3-point bending test'
rm -rf post/*
$liggghts_prog -in in.liggghts
cd ..

cd shear_box
echo 'Running shear box test'
rm -rf post/*
rm -rf restarts/*
echo 'Inserting particles'
$liggghts_prog -in in.insert_particles
echo 'Letting particles settle'
mpirun -n $num_cores $liggghts_prog -in in.settle_particles
echo 'Shearing System'
mpirun -n $num_cores $liggghts_prog -in in.shear
cd ..

cd uniaxial_compression
echo 'Running uniaxial_compression test'
rm -rf post/*
rm -rf restarts/*
echo 'Inserting particles'
$liggghts_prog -in in.write_restart
echo 'Letting particles settle'
mpirun -n $num_cores $liggghts_prog -in in.read_restart
echo 'Compressing System'
mpirun -n $num_cores $liggghts_prog -in in.compression
cd ..
