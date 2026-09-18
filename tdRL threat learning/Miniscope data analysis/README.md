# tdRL Analysis

This repository contains the MATLAB code used for the analyses reported in:

A prefrontal microcircuit computes trial-wise discount of reinforcement during threat learning 

Accepted for publication. The DOI will be added upon publication.

## Overview

The repository contains MATLAB scripts and functions for the analyses presented in the tdRL paper.


## Requirements

The analyses were performed using:

- MATLAB R2021b
- Statistics and Machine Learning Toolbox
- Optimization Toolbox

## Example data

Example data are provided in:

example.mat
run_example.m

The example dataset contains data from one example mouse and is intended to demonstrate the required MATLAB data structures.

## Input data format

### `traceeach_zs`
1 × nMice cell array

Each cell contains:
nNeurons × nFrames

Each row corresponds to one neuron and each column corresponds to one imaging frame.
`traceeach_zs` contains z-scored calcium activity for each mouse.

### `freezing_each`
1 × nMice cell array

Each cell contains:
nFrames × 1

Freezing trace for each mouse obtained using ANY-maze.

### `clustidx_each`
1 × nMice cell array

Each cell contains:
nNeurons × 1

Each value indicates the cluster assignment of the corresponding neuron.

Note: Clusters 1 and 4, based on the original cluster numbering from hierarchical clustering, are referred to as clusters 1 and 2, respectively, in the paper.

### `traceeach_S`
1 × nMice cell array

Each cell contains:
nNeurons × nFrames

Each row corresponds to one neuron and each column corresponds to one imaging frame.
`traceeach_S` contains deconvolved calcium events for each mouse.

## Running the example analysis

see `run_example.m`

The variables in `example.mat` are already organized in the format required by the analysis functions. For additional analysis, data should be organized using the same variable structures and dimensions as described.

## Setup

Download or clone this repository and add the repository and its subfolders to the MATLAB path:

addpath(genpath(pwd))

## Notes

Clusters 1 and 4, based on the original cluster numbering from hierarchical clustering, are referred to as clusters 1 and 2, respectively, in the paper.

For DRL related functions, please refer to: https://github.com/wgliee/tdRL-threat-learning.git

## Contact

For questions regarding the analysis code or data, please contact:

Mingyang Wei 
bluette0@sjtu.edu.cn