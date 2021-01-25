# LIGGGHTS Flexible Fibers

[Home](Home)

[Commands](commands)

[Tutorials](tutorial_main_page)

## How to install LIGGGHTS

Multiple steps need to be carried out to install LIGGGHTS(R).

How LIGGGHTS is installed will be greatly effected by your OS.

Please find the os-specific install instructions below.

* [Linux](how_to_install_linux)
* [Windows](how_to_install_windows)
* [Spack Installer](how_to_install_spack)

## Packages

Additional functionality can be added to the default install of liggghts by the
use of packages.

## Additional Flags

Additional flags can be used to look at additional debug information or to
change how certain models work. These flags should be added to your user
created Makefile (Makefile.user) located in src/MAKE.

There is a location for additional pre-processor flags located near the end of the user generated Makefile, under the section *LIGGGHTS pre-processor flags*.

* -DDO_SMOOTH_FIBER_FIX

    This flag adds a fix for how the bond stiffness values are calculated when
    dealing with a smooth fiber (multiple spheres overlap). This does alter the
    value of the Young's modulus and shear modulus needed to reproduce pass
    simulations. As a starting point, it is recommened to use half of the old
    modulus (If you used 2.0 GPa for the bond Young's modulus, you should now
    use 1.0 GPa).

* -DLIGGGHTS_BOND_DEBUG

    This flag turns on the debuging information associated with the bond model.

## Viewing LIGGGHTS Results

ParaView is suggested to be used when view data from LIGGGHTS (along with a few
extra plugins). For most users, please download version 5.4.1 on your native OS.

[Windows-no-mpi](https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.4&type=binary&os=Windows&downloadFile=ParaView-5.4.1-Qt5-OpenGL2-Windows-64bit.exe)

[Windows-with-mpi](https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.4&type=binary&os=Windows&downloadFile=ParaView-5.4.1-Qt5-OpenGL2-MPI-Windows-64bit.exe)
Requires [MS-MPI](https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi)

[Linux](https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.4&type=binary&os=Linux&downloadFile=ParaView-5.4.1-Qt5-OpenGL2-MPI-Linux-64bit.tar.gz)

The additional plugins that are needed can be found from the following git
repository (Please Download):

[ParaView Reader for LIGGGHTS](https://github.com/schrummy14/ParaView_Reader_for_LIGGGHTS)

Here you will find the pre-compiled plugins. The main plugin is the
'liggghts_reader'. This can be added to your ParaView by selecting:

Tools -> Manage Plugins -> Load New

This repository also holds a python macro that is used to set settings in
ParaView that make it easier to view your data. This macro is called update.py
and can be added to ParaView by:

Macros -> Add new macro
