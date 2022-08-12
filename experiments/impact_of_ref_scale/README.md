# Impact of reference scale on sourmash performance
In this experiment, we will compare the performance of sourmash on different reference scales.
There are two different scale parameters: one for the reference set, and one for the query set.
As scale increases, the data is compressed, so we expect accuracy to decrease, but speed to increase.

# Experimental setup
We will use the following parameters for all simulations:
1. 500 reference genomes
2. 1M reads for each simulation
3. 1000 genes for each simulation (`numGenes` in the `workflow.sh` script)
3. 150bp read length
4. No simulation errors

For sourmash, we will fix the following settings:
1. k=7
2. thresholdBP=100 

When analyzing the results, we will only compare filtering out 0 genes from the simulation, or else filtering out those genes that have fewer than 10 reads mapping to them.

# Execution files
The `workflow.sh` script will execute the experiment. `PlotResults.py` will visualize the results.

