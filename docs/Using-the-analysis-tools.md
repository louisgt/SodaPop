# Results and analysis

The package includes a small bioinformatics pipeline to extract, analyze and plot your data. Scripts are written in Bash and R\* and can be used separately from the main program. Those familiar with scripting languages can also use the scripts as stepping stones for further analyses.

By default, SodaPop will always output population snapshots following the parameters input by the user. This ensures you don’t lose any simulation data. To toggle on automatic analysis, use the –a [--analyze] Boolean flag with your command (see [Command line flags](command-line-flags.md)).

### \* As a requirement to make plots, the R programming language needs to be installed on your machine. You can download the [latest version of R here](https://cran.r-project.org/). The script will automatically install the required R packages for you at runtime.

## Scripts

*barcodes.sh* : this is a Bash script that uses Unix utilities such as awk, grep, uniq, join. It performs the following operations:

- a.	Convert binary snapshots to text
- b.	Extract and sort the barcodes for each time point
- c.	Compute the mean population fitness for each time point
- d.	Count the number of times a given barcode occurs at each time point
- e.	Identify the time point corresponding to the fixation of a single barcode
- f.	Combine the time points in a dataframe
- g.	Call polyclonal_structure.R with parameters

*polyclonal_structure.R* : this is a R script that imports the time series and the fitness table to make three distinct plots.

Both of these scripts can also be launched from the command-line by running barcodes.sh along with the parameters: [name of the simulation directory] [number of generations] [population size] [step (dt)] [0 if long format, 1 if short format]. As an example, the following would be the command to analyze a long format simulation with N = 10000, M = 10000 and a time step of 25 generations:

```bash
./barcodes.sh test_sim 10000 10000 25 0
```

Below you will find examples of these plots and a brief explanation.

