# SodaPop-1.0

## Overview

SodaPop is a forward-time simulator of large haploid populations aimed at studying structure, dynamics and the distribution of fitness effects without prior assumptions on the landscape. The program integrates biochemical and biophysical properties in a cell-based, object-oriented framework and provides an efficient, open-source toolkit for studying large-scale molecular evolution. SodaPop is designed with large-scale simulations in mind, making it suitable for the investigation of evolutionary dynamics in the context of antibiotic resistance, viral evolution and cancer.

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

This will use the provided makefile to build the binaries sodapop, sodasnap and sodasumm. If you wish to install the program to your computer, run the command

>
```bash
make install
```

By default, the three components above will be added to /usr/local/bin. You can change this in the makefile by editing the content of the $INSTALLDIR variable.

## API Reference

Coming soon.

## Tests

Coming soon.

## Troubleshooting

1. I get the following error when I try to compile using make: 
```
error: unrecognized command line option "-std=c++11"
```

Your gcc/g++ compiler is probably outdated. Get a newer version [here](https://gcc.gnu.org/). Anything from gcc 4.7 onwards will work.

2. I am a Windows user. How can I build and use this software?

Download and install [Cygwin](https://www.cygwin.com/) on your computer. It will allow you to build SodaPop from source and use it as you would on a Unix system.
