cd run
LIGGGHTS_TO_USE=$1
PYTHON_TO_USE=$2
$LIGGGHTS_TO_USE -in in.liggghts
$PYTHON_TO_USE results.py
rm -rf post
rm beam.csv
rm log.liggghts
cd ..
