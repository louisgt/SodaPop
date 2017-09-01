# SodaPop-1.0

## Overview

SodaPop is a forward-time simulator of large haploid populations aimed at studying structure, dynamics and the distribution of fitness effects without prior assumptions on the landscape. The program integrates biochemical and biophysical properties in a cell-based, object-oriented framework and provides an efficient, open-source toolkit for studying large-scale molecular evolution. SodaPop is designed with large-scale simulations in mind, making it suitable for the investigation of evolutionary dynamics in the context of antibiotic resistance, viral evolution and cancer.

### Visit the [SodaPop website](https://louisgt.github.io/SodaPop/) for online documentation.


## Table of Contents

[Troubleshooting](#troubleshooting)   
[Contributing](#contributing)  
[License](#license)

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
