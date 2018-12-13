# Fused Laplacian Sparse Group Lasso (FL-SGL) #

This repository contains a MATLAB implementation of the FL-SGL algorithm proposed in the paper [Modeling Alzheimer’s Disease Progression with Fused Laplacian Sparse Group Lasso](https://dl.acm.org/citation.cfm?id=3230668).

## Overview ##

The multi-task learning with fused Laplacian sparse group lasso model can model the tasks progression with a general weighted (undirected) dependency graphs among the tasks. FL-SGL encourages related tasks to have similar parameters, where the regularization depends on suitable structured sparsity based on the graph Laplacian of the task dependency matrix. FL-SGL adopts weighted task dependency graphs based on a Gaussian kernel over the time steps, which yields a fully connected graph with decaying weights. However, any task dependency graph can be used in FL-SGL. An efficient alternative directions method of multipliers based optimization algorithm is derived to solve FL-SGL formulation.

This code has been tested only in MATLAB in both Linux and Mac.

## How to run? ##

We created the file `test.m` to show how to run FL-SGL code. FL-SGL optimization considers two variants, respectively, based on multi-block ADMM and traditional two-block ADMM.

## Structure of the input data files ##

In order to run the code the input data files containing the training and test data must follow a specific format. The `FLADMM_TB()` function (two-block ADMM) and `FLADMM_MB()` function (multi-block ADMM), which are the core algorithms of FL-SGL, receive two arrays of cells, **X** (covariate matrices) and **Y** (response matrices), with both having the length equal to the number of time points *T*. In the **X** array, each cell is a matrix *n_t x d* corresponding to the input data of a specific time point. The same happens for the array **Y**, except that now each cell is a vector of size *n_t*. Note that, the number of samples at different time points can be different.


## How to cite it? ##

If you like it and want to cite it in your papers, you can use the following:

```
#!latex

@article{liu2018modeling,
  title={Modeling Alzheimer’s Disease Progression with Fused Laplacian Sparse Group Lasso},
  author={Liu, Xiaoli and Cao, Peng and Gon{\c{c}}alves, Andr{\'e} R and Zhao, Dazhe and Banerjee, Arindam},
  journal={ACM Transactions on Knowledge Discovery from Data (TKDD)},
  volume={12},
  number={6},
  pages={65},
  year={2018},
  publisher={ACM}
}
```

## Have a question? ##

If you found any bug or have a question, don't hesitate to contact me:

[Xiaoli Liu]
email: `neuxiaoliliu -at- gmail -dot- com`

