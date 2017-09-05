# Algorithm and implementation

![algo](https://user-images.githubusercontent.com/29554043/28376750-ecf38e40-6c78-11e7-92ec-3365d1dd9043.png)
**Illustration of SodaPop's core algorithm.**

## Data structures and complexity

The simulation algorithm is adapted from the Wright-Fisher model with selection. Generations are discrete time intervals in which all N parent cells within the population give birth to a certain number of daughter cells. The number of offspring k is drawn from a binomial distribution with N trials and mean w, which is the fitness of the parent cell over the sum of all cell fitnesses. The offspring go on to become the parents for the next generation.

Cell fitness is determined by the fitness function specified by the user. Among those are input-specific functions such as metabolic flux and cytotoxicity of misfolded proteins, which are dependent on protein stability (∆G). The distribution of fitness effects (DFE) can be input explicitly by the user in the form of a ∆∆G matrix or of deep mutational scanning (DMS) data. Fitness effects can also be drawn from a Gaussian distribution with a given mean and SD.

The population is implemented as a vector (*std::vector* container) of cells. This is because vectors store elements contiguously in physical memory, making elements accessible sequentially through iterators and pointer arithmetic, and providing spatial locality of reference. The most common operation is appending a cell at the end of the vector, which is done in amortized constant time. This is linear if we append k cells at once. Random iterator access is also done in constant time, making the random rescaling of the population an efficient operation. Vectors also confer efficient memory handling as memory can be reserved ahead of time, provided a good estimate of capacity is known. Because we have fixed size populations, this is a known variable, preventing costly vector doubling reallocations when the container is at full capacity. Cells also use vectors to implement genomes. In this implementation, the genomic structure is static. There are thus no operations other than iteration applied on cell vectors.

Whenever it is possible, we use safe C++ standard library operators such as *std::shuffle* and *std::swap*. The first operator is used when the new generation exceeds the intended size (as a result of binomial drawings). It performs a fast random shuffle of the population vector. This is followed by a standard resizing of the vector whose complexity is linear on the number of elements erased. The second operator is used to assign the new vector of cells to the population vector, overwriting the previous generation. Because cells are deep data structures with a lot of information, overwriting the population one by one would not only be inefficient, but potentially pose the risk of losing information by shallow copies. The *swap* operator uses *move* semantics to swap the two cell vectors efficiently. 

The program reads and writes binary files because they are smaller in size. Since the memory overhead for I/O is significantly higher than that for type conversion, the computational cost incurred by converting types to binary is negligible. Moreover, while using binary over raw text might not be significant for small-scale simulations, it can trim down the size of snapshots in larger runs where each population can reach several hundred megabytes in size.  

## Pseudo-random number generation

We are working with the assumption that arising mutations are uniformly distributed along the genome, and a single simulation run can easily produce millions of mutations. Hence, we need a reliable pseudo-random number generator (PRNG) that lacks bias. Additionally, we want our PRNG to be relatively fast as it is being used heavily in the algorithm. Finally, the PRNG should be lightweight and compatible with C/C++. Following these criteria, we opted to use the PCG family of random generators by Melissa E. O’Neill, as a substitute to the widely used Mersenne Twister provided in C++11. The reason is that the latter has a large state space and is known to be uneven in its output, while PCG is compact, relatively fast and uniform in its output, which makes it a good candidate for simulation.

## Go back to the [home page.](index.md)
