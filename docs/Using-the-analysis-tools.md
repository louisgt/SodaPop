# Results and analysis

The package includes a small bioinformatics pipeline to extract, analyze and plot your data. Scripts are written in Bash and R\* and can be used separately from the main program. Those familiar with scripting languages can also use the scripts as stepping stones for further analyses.

By default, SodaPop will always output population snapshots following the parameters input by the user. This ensures you don’t lose any simulation data. To toggle on automatic analysis, use the –a [--analyze] Boolean flag with your command (see [Command line flags](command-line-flags.md)).

### \* As a requirement to make plots, the R programming language needs to be installed on your machine. You can download the [latest version of R here](https://cran.r-project.org/). The script will automatically install the required R packages for you at runtime.

## Scripts

***barcodes.sh*** : this is a Bash script that uses Unix utilities such as awk, grep, uniq, join. It performs the following operations:

- a.	Convert binary snapshots to text
- b.	Extract and sort the barcodes for each time point
- c.	Compute the mean population fitness for each time point
- d.	Count the number of times a given barcode occurs at each time point
- e.	Identify the time point corresponding to the fixation of a single barcode
- f.	Combine the time points in a dataframe
- g.	Call polyclonal_structure.R with parameters

***polyclonal_structure.R*** : this is a R script that imports the time series and the fitness table to make three distinct plots.

Both of these scripts can also be launched from the command-line by running barcodes.sh along with the parameters: [name of the simulation directory] [number of generations] [population size] [step (dt)] [0 if long format, 1 if short format]. As an example, the following would be the command to analyze a long format simulation with N = 10000, M = 10000 and a time step of 25 generations:

```bash
  ./barcodes.sh test_sim 10000 10000 25 0
```

Below you will find examples of these plots and a brief explanation.

## Mean fitness plot

The first plot is titled *fitness.tiff* and represents the curve of the mean population fitness during the course of the simulation.

![fitness example](https://user-images.githubusercontent.com/29554043/29976715-ed9fea34-8f08-11e7-82be-d8800e4ec475.png)

In the example above, we see a steep drop in fitness at the start. This is a classic case of Muller’s ratchet: under a high mutation rate, the population sees its mean fitness decrease as the fitter individuals are hit by marginally deleterious mutations. Over time, marginally beneficial mutations confer a selective advantage that drives the mean fitness higher.

## Clonal structure plot

The second plot is titled clonal_structure.tiff and represents the fraction of the population held by each segregating lineage, that is, cells sharing a common barcode. 

The vertical axis represents the density of each lineage. The sum of all lineages is equal to the size of the population (in the case above, 10,000). Attentive observers may notice the slope in the top left area of the plot. This is a result of discarding all lineages that are lost immediately after the first time point, the reason being that as the ratio of unique barcodes to the population size approaches 1, a majority of barcodes will be naturally lost to stochastic drift. Removing these lineages speeds up plotting and makes for a sharper image.

![fitness example](https://user-images.githubusercontent.com/29554043/29976704-e4433676-8f08-11e7-9421-a02f6dad4e98.png)

By default, this plot will show lineages up to the moment of fixation (if it exists), as the subsequent time points would only show a static image.

## Clonal trajectories plot

The third and last plot is titled clonal_trajectories.tiff. It is an alternative visualization of the previous plot.

![fitness example](https://user-images.githubusercontent.com/29554043/29976708-e9a47558-8f08-11e7-9069-9195e4accc87.png)

The vertical axis represents the total count of each lineage as it segregates. Again, this plot will show lineages up to the moment of fixation.

## Mutation log file

The file MUTATION_LOG contains the information on arising mutations in a simulation run. It can be toggled on with the –e [--track-events] flag.

The tab-separated file lists mutations line-by-line:

```
AGACTCAAGTGTGAC 6       K       272     E       0.001078        1
TGCAAGCAAACGGGC 5       D       160     H       -0.363828       9
CGAGACTGTGGGAGT 4       T       48      R       -0.198414       21
...
```

The columns respectively correspond to the barcode, the gene ID, the prior amino acid, its position in the sequence, the resulting mutation, its selection coefficient and the generation at which it occurred. 


# Move on to [command-line flags](command-line-flags.md).

## or go back to [home page.](index.md)
