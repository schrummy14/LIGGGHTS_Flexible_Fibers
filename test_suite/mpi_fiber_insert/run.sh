cd run
LIGGGHTS_TO_USE=$1
mpirun -n 4 $LIGGGHTS_TO_USE -in in.liggghts
rm -rf post
rm log.liggghts
cd ..
