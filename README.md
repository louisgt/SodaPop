# SodaPop-1.0

## Developed by the http://www.serohijoslab.org/

## Visit the [SodaPop website](https://louisgt.github.io/SodaPop/) for online documentation or [download the SodaPop user manual](https://github.com/louisgt/SodaPop/files/2036999/manual.pdf).

## Overview

SodaPop is a forward-time simulator of large asexual populations aimed at studying population structure, dynamics and the distribution of fitness effects without prior assumptions on the landscape. The program integrates biochemical and biophysical properties in a cell-based, object-oriented framework and provides an efficient, open-source toolkit for studying multi-scale scenarios in molecular evolution.

![header](https://user-images.githubusercontent.com/29554043/56760304-e6537680-6768-11e9-9f40-21010e5b3c93.png)

## Installation

To decompress and extract the contents of the downloaded repository, open a command-line terminal window, change into the directory where the download is located on your computer and run the following command

```bash
  tar –zxvf [zip file]
```

Before you proceed with the installation, you may wish to move the extracted files to a folder of your choosing. To compile the SodaPop locally, navigate to the program folder and run

```bash
  make
```

This will use the makefile to build the binaries **sodapop**, **sodasnap** and **sodasumm**. To install the program to your computer, run the command


```bash
  make install
```

By default, the three components above will be added to /usr/local/bin. You can change this in the makefile by editing the content of the $INSTALLDIR variable. Likewise, any other parameter in the makefile can easily be modified.

## Troubleshooting

1. I get the following error when I try to compile using make: 
```
error: unrecognized command line option "-std=c++11"
```

Make sure your gcc/g++ compiler is up-to-date. Get a newer version [here](https://gcc.gnu.org/). Anything from gcc 4.7 onwards will work.

2. I am a Windows user. How can I build and use this software?

Download and install [Cygwin](https://www.cygwin.com/) on your computer. This will allow you to build SodaPop from source and use it as you would on a Unix system.


## Issues and bugs

SodaPop has been tested extensively on a variety of machines. We have taken great care to test diverse scenarios and to correct any bugs we could find. With that being said, if you stumble upon bugs, compilation or execution issues, or if you have any questions on the software, please do share those with us, either by [creating an issue](https://help.github.com/en/articles/creating-an-issue) on the project page or by email at louis.gauthier AT umontreal.ca .

## Contributing

SodaPop is a work in progress. It can certainly be expanded by implementing new features, and it can also be further optimized for performance and memory. If you have ideas for additional features, or if you are interested in contributing to the software, please contact me at louis.gauthier AT umontreal.ca

## Citation

*If you use SodaPop in your research, kindly cite the following:*

Gauthier L, Di Franco R, Serohijos AWR. SodaPop: A Forward Simulation Suite for the Evolutionary Dynamics of Asexual Populations on Protein Fitness Landscapes. Bioinformatics. 2019 Mar 13. pii: btz175. doi:10.1093/bioinformatics/btz175.

```
@article{Gauthier2019,
    author = {Gauthier, Louis and Serohijos, Adrian W R and Di Franco, Rémicia},
    title = "{SodaPop: A Forward Simulation Suite for the Evolutionary Dynamics of Asexual Populations on Protein Fitness Landscapes}",
    year = {2019},
    doi = {10.1093/bioinformatics/btz175}
}
```


## License

Copyright (C) 2019 Louis Gauthier

SodaPop is free software. You can redistribute it or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

SodaPop is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with SodaPop.  If not, see <http://www.gnu.org/licenses/>.
