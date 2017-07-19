To create a binary population snapshot ‘from scratch’, you can use the program called sodasumm. It requires a population description file, which consists of the count of each cell type to be included in the population. The program can generate either a full population snapshot or the information for a single cell. The latter may be preferable in the case of simulations involving very large populations. Starting with a monoclonal population effectively speeds up the initialization process, as the vector is populated with copies of the same cell. In the former case, cells must be added one at a time, which implies a cost that scales linearly with the number of cells. Nonetheless, it is possible for the user to choose either option.

> $ sodasumm <population summary> [ 0-full | 1-single cell ]

The population summary file simply lists the paths to the cell files and their respective counts.

![Pop_dat](https://user-images.githubusercontent.com/29554043/28380453-61b90df2-6c85-11e7-881a-e964445f068c.png)

