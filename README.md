# SodaPop-1.0

## Developed by the http://www.serohijoslab.org/

## Visit the [SodaPop website](https://louisgt.github.io/SodaPop/) for online documentation or [download the SodaPop user manual](https://github.com/louisgt/SodaPop/files/1271058/MANUAL.pdf).

## Overview

SodaPop is a forward-time simulator of large asexual populations aimed at studying population structure, dynamics and the distribution of fitness effects without prior assumptions on the landscape. The program integrates biochemical and biophysical properties in a cell-based, object-oriented framework and provides an efficient, open-source toolkit for studying large-scale scenarios in molecular evolution.

![header](https://user-images.githubusercontent.com/29554043/32801437-0cb5cc14-c94b-11e7-8b22-5687ff245afc.png)

## Troubleshooting

1. I get the following error when I try to compile using make: 
```
error: unrecognized command line option "-std=c++11"
```

Make sure your gcc/g++ compiler is up-to-date. Get a newer version [here](https://gcc.gnu.org/). Anything from gcc 4.7 onwards will work.

2. I am a Windows user. How can I build and use this software?

Download and install [Cygwin](https://www.cygwin.com/) on your computer. This will allow you to build SodaPop from source and use it as you would on a Unix system.


## Contributing

SodaPop is a work in progress. It can certainly be expanded by implementing new features, and it can also be further optimized for performance and memory. If you have ideas for additional features, or if you are interested in contributing to the software, please contact me at louis.gauthier AT umontreal.ca


## License

Copyright (C) 2018 Louis Gauthier

SodaPop is free software. You can redistribute it or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

SodaPop is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with SodaPop.  If not, see <http://www.gnu.org/licenses/>.
