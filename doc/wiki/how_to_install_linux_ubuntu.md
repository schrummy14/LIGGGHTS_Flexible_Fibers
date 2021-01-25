# LIGGGHTS Flexible Fibers

[Home](Home)

[Commands](commands)

[Tutorials](tutorial_main_page)

[How To Install](how_to_install)

## How to Install LIGGGHTS on Ubuntu

* Install Packages

```text
sudo apt-get install build-essential libopenmpi-dev libvtk7-dev
```

* Clone the repository

```text
git clone https://github.com/schrummy14/LIGGGHTS_Flexible_Fibers
```

* Make the package

```text
cd src
make -j4 auto
```

* Create a symbolic link

```text
sudo ln -s ~/where/you/installed/liggghts/src/lmp_auto /usr/local/bin/liggghts_bonds
```

* Run the test suite (The "bond break" test may fail and this is ok as it is
experimental)

```text
cd ..
cd test_suite
python3 testAll.py
```
