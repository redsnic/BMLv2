# BMLv2

Unofficial updated version of BML:Bayesan Mutation Landscape by Navodit Misra et al., includes a newly made GUI (made with the [Qt library](https://www.qt.io/)), a MAF preprocessor and multithreading support. It is compatible with both Windows and Linux.

# How to get Qt libraries

It is possible to download Qt directly from [here](https://www1.qt.io/download-open-source). Make sure that your version is compatible with *webenginewidgets* and *webengine*, tested version are 5.9.1 (gcc) on Linux and 5.8.0 (msvc2013) on Windows.

# How to Compile on Linux

To compile this programm you need to download the Source folder on your computer, move yourself inside it and run:

`qmake` or `./PathToYourQtInstallation/bin/qmake`

and then

`make`.

You can now execute the program: 

`./BML_GUI`

> I haven't uploaded the Windows .pro file needed by Qt yet, coming soon.

# Changelog
## v1.0
First publicated version by Navodit Misra et al., original sources can be downloaded from [here](http://bml.molgen.mpg.de/)
## v2.0
This is the last version. Was prepared by Nicolò Rossi (nicolo@rossi.im) and consists in these additions:
* a MAF file preprocessor
* a complete Graphical User Interface
* multithreading options

The code was also partially revisited by adding various comments, changing variables names with more informative ones and solving some minor bugs.
