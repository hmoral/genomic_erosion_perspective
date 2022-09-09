# -*- coding: utf-8 -*-
"""
Script burnin from msprime
For simulations in van Oosterhout et al. 2022 preprint: Genomic erosion in the assessment of species extinction risk and recovery potential
written by Hernán E. Morales (hern.moral@gmail.com) 
"""
import msprime, pyslim, tskit
import numpy as np
import argparse
#import allel

print(msprime.__version__)
print(pyslim.__version__)
print(tskit.__version__)
#print(allel.__version__)


# =============================================================================
# outPath="./0_github/"
# outName="out_uncondLoad_msprimeBurn_exomeONLY"
# geneLength=3400
# g=1000
# U=0.1
# Ne=1000
# seed=12345
# chrNum=10;
# 
# =============================================================================

#### Get arguments
parser = argparse.ArgumentParser()
# tree name
# Ne for recapitation
parser.add_argument('--Ne',nargs='+',type=str)
# recombination rates file
# U_in inoput mutation rate exome
parser.add_argument('--U',nargs='+',type=float)
# exome gene len
parser.add_argument('--geneLength',nargs='+',type=str)
# exome gene count
parser.add_argument('--g',nargs='+',type=str)
# number of chromosomes
parser.add_argument('--chrNum',nargs='+',type=int)
# output name
parser.add_argument('--outName',nargs='+',type=str)
# output path
parser.add_argument('--outPath',nargs='+',type=str)
# seed
parser.add_argument('--seed',nargs='+',type=int)
# recombination rate
parser.add_argument('--rho',nargs='+',type=float)


args = parser.parse_args()
Ne=float(args.Ne[0])
U=float(args.U[0])
geneLength=int(args.geneLength[0])
g=int(args.g[0])
chrNum=int(args.chrNum[0])
outName=str(args.outName[0])
outPath=str(args.outPath[0])
seed=int(args.seed[0])
rho=float(args.rho[0])

Glen=float(geneLength*g)
mu=U/(2*Glen)
mu_exome=mu*(1/3.31)

print("Exome len:"+str(Glen/1e6)+"Mb")
print("U:"+str(U))

print("mutation rates:")
print("mutation rate total:"+str(mu))
print("mutation exome:"+str(mu_exome))
print("recombination rate:")
print(str(rho))

# =============================================================================
# mu_rates=[]
# mu_ends=[]
# for i in range(1,g):
#     mu_rates.append(mu)
#     mu_ends.append((i*geneLength)+(i-2))
# 
# =============================================================================

gene_nums=list(np.tile(g/chrNum, chrNum))
rates=[];
rho=1e-4;
for i in range(1,len(gene_nums)):
    print(i)
    rates.append(0)
    rates.extend(list(np.tile([rho, 0], int(gene_nums[i-1]-1))))
    rates.append(0.5)
rates.append(0)
rates.extend(list(np.tile([rho, 0], int(gene_nums[len(gene_nums)-1]-1))))


ends=[];
for i in range(1,g+1):
    ends.append((i*geneLength)+(i-2))
    ends.append((i*geneLength)+(i-1))
ends=ends[0:(len(ends)-1)]
ends.insert(0,0)
ends[-1]=ends[-1]+1


totalLen=ends[len(ends)-1]

recomb_map = msprime.RateMap(
  position = ends,
  rate = rates)

demog_model = msprime.Demography()
demog_model.add_population(initial_size=Ne)
ots = msprime.sim_ancestry(
        samples=Ne,
        demography=demog_model,
        random_seed=seed,
        recombination_rate=recomb_map)

ots = pyslim.annotate_defaults(ots, model_type="nonWF", slim_generation=1)


GENS=int(Ne*10000)
mut_model = msprime.SLiMMutationModel(type=1,slim_generation=GENS,next_id=1)              
ots = msprime.sim_mutations(
            ots,
            rate=mu_exome,
            model=mut_model,
            keep=True, 
            random_seed=seed)
print(f"The tree sequence now has {ots.num_mutations} mutations, at "
      f"{ots.num_sites} distinct sites.")


#FILEname=str(outPath)+"/"+str(outName)+"_g"+str(g)+"_geneLength"+str(geneLength)+"_muChrom"+str(muChrom_in)+"_Ne"+str(Ne)+"_U"+str(U_in)+"_seed"+str(seed)+"_burn.tree" 
#print('Output '+FILEname)
#ots.dump(FILEname)

genotypes=np.array(ots.genotype_matrix())
mut_count=np.sum(genotypes, axis = 1)
pi_diver=float(ots.diversity())

print("pi="+str(pi_diver))


FILEname=str(outPath)+"/"+str(outName)+"_g"+str(g)+"_geneLength"+str(geneLength)+"_totalLen"+str(totalLen)+"_Ne"+str(Ne)+"_U"+str(U)+"_muNeutralExome"+str(mu_exome)+"_seed"+str(seed)+"_burn.tree" 
print('Output '+FILEname)
ots.dump(FILEname)
