# LIGGGHTS Flexible Fibers

[Home](Home)

[Commands](commands)

[Tutorials](tutorial_main_page)

## How to install LIGGGHTS

---

## Linux

---

### Ubuntu

1) Install Packages

```text
sudo apt install ...
```

2) Clone the repository

```text
git clone https://github.com/schrummy14/LIGGGHTS_Flexible_Fibers
```

3) Make the package

```text
cd src
make -j4 auto
```

4) Create a symbolic link

```text
sudo ln -s ~/where/you/installed/liggghts/src/lmp_auto /usr/local/bin/liggghts_bonds
```

5) Run the test suite

```text
cd ..
cd test_suite
python3 testAll.py
```

---

## Windows

---

### WSL

This should be the same as the distribution inside the WSL

### Native

This is can be fun...
