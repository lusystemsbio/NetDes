# NetDes
Network inference and optimization using Dynamical equation simulations

NetDes is a computational method for modeling gene regulatory networks (GRNs) using time-series scRNA-seq data. It infers pseudotime, identifies gene expression trajectories, clusters genes, and detects key transcription factors (TFs). A GRN is constructed and refined using ordinary differential equations (ODEs), enabling dynamic simulations and analysis of GRN-driven cell state transitions.

## Installation
To install NetDes, Python version 3.9 or greater is required.

### Install from PyPi (recommended)
To install the most recent release, run

`pip install NetDes`

### Install with github
* Git clone the [NetDes repository](https://github.com/lusystemsbio/NetDes), cd to the `NetDes` directory, and run

`pip install .`

## Tutorials

### Prepare

[Data Pre-processing, Gene Expression Smoothing and Clustering ](tutorials/R_dataprocess/1_Trajectores_and_clusters.html): This tutorial walks through the process of inputting scRNA-seq data, processing it, and obtaining smooth gene expression trajectories and clustering them.

[TFs Identification](tutorials/R_dataprocess/2_TFs_identify.html): This tutorial demonstrates how to use Fisher's exact test to identify core transcription factors.

### Network inferring and downstream analysis
[Initial GRN Building](tutorials/R_dataprocess/3_InitialGRN.html): Learn how to load databases like Rcistarge, TRRUST, and NetAct to construct the initial gene regulatory network.


[GRN Inferring and Simulation](tutorials/tutorial.html): This tutorial covers how to use pseudotime information, smoothed gene expression trajectories, and the initial GRN for inference and simulation. It includes a practical example using synthetic data, which you can generate with [this script.](tutorials/datasimulation.html)

[Network Evaluation](tutorials/R_dataprocess/4_GRN_evaluation.html): After building the optimized GRN, this tutorial explains how to evaluate its performance through simulations.

[Coarse-graining the Network](tutorials/R_dataprocess/5_Coarse_graining.html): This tutorial walks through the process of coarse-graining the optimized GRN into a smaller circuit.


