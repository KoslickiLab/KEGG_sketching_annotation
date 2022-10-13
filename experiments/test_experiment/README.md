This is a trial experiment to make sure that things are working.

# Experiment setup
We will use the following parameters for all simulations:
1. 30 reference genomes, genes from 25 will be present in ref db, genes from the
rest 5 will not be present in the ref db. The number 25 ensures that there are around 100k genes in the ref db, so that sourmash gather does not take forever
1. 1M reads for each simulation
1. 10000 genes for each simulation
1. 150bp read length
1. No simulation errors

For sourmash, we will fix the following settings:
1. k=7
1. thresholdBP=100
