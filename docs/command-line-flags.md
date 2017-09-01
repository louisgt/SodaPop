
# Command-line flags

Command-line flags are modifiers that toggle features, change behaviors and provide great flexibility on the simulation parameters. A description of the flags can always be displayed in the command-line itself by typing 

```
sodapop --help
```

Flags can also be added or modified in the evolve.cpp source file, using the [TCLAP library API](http://tclap.sourceforge.net/manual.html).

```
m [maxgen]: this is an integer specifying the length of the simulation
```

```
n [size]: this is an integer defining the target size of the population following a reproduction step
```

```
t [dt]: this is an integer specifying the interval at which a population snapshot is to be produced. The minimum value is 1. There is no maximum. If the interval exceeds the number of generations, only one snapshot will be output, prior to the evolutionary process.
```

```
o [prefix]: this is the prefix that will be used for output. In case the specified directory does not exist, it will be created. Otherwise, it will be overwritten. This parameter is not required but nevertheless strongly recommended as the default setting will create a directory named ‘sim/’ and possibly overwrite information of previous runs. You should always specify a new output if you want to distinguish between multiple runs.
```

```
g [gene-list]: this is the file that lists all the genes (and their respective file name) pertaining to the simulation.
```

```
p [pop-desc]: this is the binary snapshot used to initialize the population in the simulation.
```

```
l [gene-lib]: this is the path to the folder containing the gene files. This option is not necessary if you put your genes in the existing folder (/files/genes/). You should specify the path if you use a different folder.
```

```
i [input]: this is the file that defines the fitness landscape for the specified genes. It contains a matrix of values for each different gene. It may exceed the actual number of genes used in simulations as long as the gene numbers are consistent with the .gene files. The format is the same for all input types, the default type being deep mutational scanning (DMS) selection coefficients. Other types are specified with the appropriate flags (see below). For a given gene matrix M, the value corresponding to the substitution of amino acid a at site k is indexed as M_a,k. For simplicity, all 20 amino acid substitutions are listed in alphabetical order (from ‘A’ to ‘Y’). The internal mappings of the genetic code are represented this way. The identity substitution (from the native amino acid to itself) must be included and given a neutral numeric value.
```

```bash
f [fitness]: the user can select from several fitness functions depending on the input type and the biophysical or biochemical property of interest. For a detailed reference of available functions, refer to section 6.
```

```
sim-type: this flag is used to specify the theoretical and/or experimental background that will guide the evolutionary process. The default type is deep mutational scanning (DMS) selection coefficient. Alternatively, the type can be set to thermodynamic stability (∆G) or another property implemented by the user.
```

```
gamma: toggling this flag will activate drawing values from a gamma distribution. The shape and scale parameters are specified with the alpha and beta flags (see below).
```

```
normal: toggling this flag will activate drawing values from a normal distribution. The mean and standard deviation are specified with the alpha and beta flags (see below).
```

```
alpha: floating-point value of the alpha parameter. If gamma distribution is used, this parameter corresponds to the shape. If normal distribution is used, this parameter corresponds to the mean.
```

```
beta: floating-point value of the beta parameter. If gamma distribution is used, this parameter corresponds to the scale. If normal distribution is used, this parameter corresponds to the standard deviation.
```

```
c [create-single]: toggling this flag will populate the initial vector with duplicates of a single cell. This can dramatically speed up initialization for large-scale populations, at the cost of clonal diversity. If your starting population is polyclonal, you should not use this flag.
```

```
a [analyze]: toggling this flag will call analysis scripts after the simulation has ended.
```

```
e [track-events]: toggling this flag will keep a log of all arising mutations during the course of the simulations.
```

```
s [short-format]: toggling this flag will save a more compact snapshot of the population to track population dynamics only. As mentioned above, normal snapshot size scales linearly with the number of cells and the size of their genome. Thus, large simulations can output an equivalently large volume of data. The information for individual genes (including explicit sequences) is not saved when using this format. If sequence information is not needed for an experiment, users may choose to toggle this option to minimize output. Unless disk space is scarce, we strongly recommend keeping the normal output format when running simulations. This will ensure you do not discard valuable information.
```

# Move on to [description of the algorithm](Description-of-the-core-algorithm.md).

## or go back to the [home page.](index.md)
