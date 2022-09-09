### The code is run in three steps:
1. extRisk_perspective_WFnonWF_msprime_uncondLoad_NeutralBurnIN_exomeONLY.py: used for generating a burnin stage; produces a mutated treeseq file using msprime (see https://tskit.dev/msprime/docs/stable/intro.html)
2. extRisk_perspective_WFnonWF_msprime_uncondLoad_tree2Full.slim: transform the tree from step 1 into slim output, the tree can be read directly into step 3 but the simulation is slower because it needs treeseq recording
3. extRisk_perspective_WFnonWF_msprime_uncondLoad_exomeONLY.slim: takes mutated output from step 2 and runs the simulation of unconditional genetic load 

### minimal example (unrealistic parameters, e.g. small population and few genes to run faster) 
```python extRisk_perspective_WFnonWF_msprime_uncondLoad_NeutralBurnIN_exomeONLY.py --Ne 1000 --U 0.1 --geneLength 1000 --g 1000 --chrNum 10 --outName out_minimal_test --outPath ./ --seed 1234 --rho 0.0001```

```slim -d "K=1000" -d "outPath='./'" -d "treeIN='out_minimal_test_g1000_geneLength1000_totalLen1000999_Ne1000.0_U0.1_muNeutralExome1.5105740181268882e-08_seed1234_burn.tree'" -d "totalLen=1000999" extRisk_perspective_WFnonWF_msprime_uncondLoad_tree2Full.slim```

```slim -d "K=1000" -d "outPath='./'" -d "outName='out_minimal_test'" -d "neutralBurn='out_minimal_test_g1000_geneLength1000_totalLen1000999_Ne1000.0_U0.1_muNeutralExome1.5105740181268882e-08_seed1234_burn_fullOut.slim'" -d "rho=0.0001" -d "seed=1234" -d "U=0.1" -d "geneLength=1000" -d "g=1000" -d "chrNum=10" extRisk_perspective_WFnonWF_msprime_uncondLoad_exomeONLY.slim```

Dendencies:
SLiM3, python and python libraries: msprime, pyslim, tskit, numpy, argparse, allel and re
