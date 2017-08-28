# SodaPop-1.0

## Overview

SodaPop is a forward-time simulator of large haploid populations aimed at studying structure, dynamics and the distribution of fitness effects without prior assumptions on the landscape. The program integrates biochemical and biophysical properties in a cell-based, object-oriented framework and provides an efficient, open-source toolkit for studying large-scale molecular evolution. SodaPop is designed with large-scale simulations in mind, making it suitable for the investigation of evolutionary dynamics in the context of antibiotic resistance, viral evolution and cancer.


## Table of Contents

[Installation](#installation)  
[Usage](#usage)  
[Flags](#flags)  
[Troubleshooting](#troubleshooting)   
[Contributing](#contributing)  
[License](#license)

<a name="installation"/>

## Installation

SodaPop is set of command-line tools written in the C++ programming language using the C++11 Standard. To decompress and extract the contents of the downloaded repository, open a command-line terminal window, change into the directory where the download is located on your computer and run the following command

>
```bash
tar –zxvf [zip file]
```

Before you proceed with the installation, you may wish to move the extracted files to a folder of your choosing. To compile the SodaPop locally, navigate to the program folder and run

>
```bash
make
```

This will use the makefile to build the binaries **sodapop**, **sodasnap** and **sodasumm**. To install the program to your computer, run the command

>
```bash
make install
```

By default, the three components above will be added to /usr/local/bin. You can change this in the makefile by editing the content of the $INSTALLDIR variable. Likewise, any other parameter in the makefile can easily be modified.

<a name="usage"/>

## Usage

For a complete walkthrough of commands, scripts and output, please refer to the [SodaPop wiki](https://github.com/louisgt/SodaPop/wiki).

![Pop. dynamics example](https://user-images.githubusercontent.com/29554043/28281174-42b56b7c-6af4-11e7-86c9-f8393c123513.png)

![Fitness trajectory](https://user-images.githubusercontent.com/29554043/28281203-573643f0-6af4-11e7-9362-212a833a056f.png)

![Clonal interference](https://user-images.githubusercontent.com/29554043/28374005-7e0a6e60-6c71-11e7-945b-546541b0c3f1.png)

<a name="flags"/>

## Flags

>
```bash
m [maxgen]: number of generations to simulate
```

>
```bash
n [size]: population size cutoff
```

>
```bash
t [dt]: time interval for population snapshots
```

>
```bash
o [prefix]: prefix to use for output
```

>
```bash
g [gene-list]: gene list file
```

>
```bash
p [pop-desc]: population description (snapshot) file
```

>
```bash
l [gene-lib]: path to gene library
```

>
```bash
i [input]: input file defining the fitness landscape
```

>
```bash
f [fitness]: fitness function to use
```

>
```bash
sim-type: defines the theoretical and experimental background for the simulation
```

>
```bash
gamma: draw values from gamma distribution
```

>
```bash
normal: draw values from normal distribution
```

>
```bash
alpha: alpha parameter of distribution
```

>
```bash
beta: beta parameter of distribution
```

>
```bash
c [create-single]: create initial population from a single cell
```

>
```bash
a [analyze]: enable automatic analysis of simulation
```

>
```bash
e [track-events]: track arising mutations
```

>
```bash
s [short-format]: use short format for population snapshots
```

<a name="troubleshooting"/>

## Troubleshooting

1. I get the following error when I try to compile using make: 
```
error: unrecognized command line option "-std=c++11"
```

Make sure your gcc/g++ compiler is up-to-date. Get a newer version [here](https://gcc.gnu.org/). Anything from gcc 4.7 onwards will work.

2. I am a Windows user. How can I build and use this software?

Download and install [Cygwin](https://www.cygwin.com/) on your computer. This will allow you to build SodaPop from source and use it as you would on a Unix system.

<a name="contributing"/>

## Contributing

SodaPop is a work in progress. It can be expanded by implementing new features. It can certainly be further optimized for performance and memory. If you have ideas for additionnal features, or if you are interested in contributing to the software, please contact Louis Gauthier at louis.gauthier@umontreal.ca

<a name="license"/>

## License

Copyright (C) 2017 Louis Gauthier

SodaPop is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

SodaPop is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with SodaPop.  If not, see <http://www.gnu.org/licenses/>.
