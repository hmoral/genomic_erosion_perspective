
## Instructions for running simulations for hard and soft selection
Both the hard and soft selection simulations are run in two steps using SLiM (URL for SLiM manual - http://benhaller.com/slim/SLiM_Manual.pdf)

Both simulations require the the generation of 'burn in' populations that at generated using "Burn_in_populations_sript.slim"
It requires deleterious mutations sampled from previously simulated stable populations which are read in using "hs_in.txt", and data used to create the exomes from the populations are read in from "chr_genes_in.txt".
The paths to each of these files needs completed within the SLiM scripts.
Paths that to the desired destination for the output of data also need to be completed.
Each simulation ran using "Burn_in_populations_sript.slim" produces a population file and a file of selected mutations.
Both of these files needs to be read into the hard and soft selection simulations in parts 2a and 2b respectively, in addition to the "chr_genes_in.txt" file.

- 2a The simulations that produce data for genetic load under hard selection are run using "Hard_selection_script.slim"
- 2b The simulations that produce data for genetic load under soft selection are run using "Soft_selection_script.slim"
These scripts also require the completion of the paths used to read files.
Paths that to the desired destination for the output of data also need to be completed.

## Command line example

```slim -s 135 -d K=5000 -d Ls=5 -d hmode=2 -d GenLoadth=1.7 -d geneLength=1000 -d g=1000 -d rho=1e-05 -d mu=1e-05 Hard_selection_script.slim```

Key:

- s           Seed number
- K            The mean of the normal distribution from which the carrying capacity will be selected each generation
- Ls           The mean of the poisson distribution from which the number of offspring from each mating couple will be selected
- hmode=2      Ensures that mutations are selected from the relevant files rather than generated independently
- GenLoadth    The threshold that the each individual needs to score over to survive hard selection (hard selection simulations only)
- geneLength   The number of base pairs of each gene
- g            The number of genes
- rho          Recombination rate
- mu           Mutation rate
