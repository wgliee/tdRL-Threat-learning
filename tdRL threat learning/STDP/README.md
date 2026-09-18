# mPFC microcircuit model

This folder contains the MATLAB code used to model the mPFC microcircuit shown in Fig. 6r-w. The model compares Hebbian spike-timing-dependent plasticity (STDP), anti-Hebbian STDP, and a no-STDP condition using leaky integrate-and-fire neurons.

## Requirements

- MATLAB R2021b
- Statistics and Machine Learning Toolbox (`poissrnd`)

No external input data or third-party packages are required.

## Files

- `simplifiedmodel_hebbian.m`: Hebbian STDP condition
- `simplifiedmodel_anti.m`: anti-Hebbian STDP condition
- `simplifiedmodel_noSTDP.m`: no-STDP condition

## Running the model

Open MATLAB, set the current folder to this directory, and run one of the following commands:

```matlab
V_hebbian = simplifiedmodel_hebbian;
V_anti = simplifiedmodel_anti;
V_noSTDP = simplifiedmodel_noSTDP;
```

Each function runs one stochastic simulation and returns `V`, a `1000 x 100` matrix of membrane potentials. Rows correspond to simulation time points. Columns correspond to neurons:

- columns 1-40: mPFC SST interneurons;
- columns 41-70: cluster 1 excitatory neurons; and
- columns 71-100: cluster 2 excitatory neurons.

Because connectivity and Poisson noise are generated randomly, use `rng` before calling a function when a repeatable example is needed:

```matlab
rng(1, 'twister');
V_hebbian = simplifiedmodel_hebbian;
```

To reproduce the averaging procedure described in the manuscript, run each condition 100 times and average the returned matrices:

```matlab
nRepetitions = 100;
V_hebbian_all = zeros(1000, 100, nRepetitions);

for repetition = 1:nRepetitions
    V_hebbian_all(:, :, repetition) = simplifiedmodel_hebbian;
end

V_hebbian_mean = mean(V_hebbian_all, 3);
```

Repeat the same procedure with `simplifiedmodel_anti` and `simplifiedmodel_noSTDP`.

## Model parameters

The code follows the parameters reported in the Methods:

- 100 neurons: 40 interneurons, 30 cluster 1 neurons, and 30 cluster 2 neurons;
- resting potential: -60 mV;
- firing threshold: -40 mV;
- firing potential: 20 mV;
- reset potential: -70 mV;
- membrane time constant (`C/g`): 10 time steps;
- interneuron connection probability: 0.2 for cluster 1 and 0.5 for cluster 2;
- excitatory-to-excitatory connection probability: 0.1;
- excitatory and inhibitory synaptic weights: 0.03 and -0.01;
- STDP window: five time steps;
- maximum STDP weight adjustment: 10%; and
- five external inputs with increasing interneuron input strength.

In the weight matrix, `1` and `-1` indicate the presence of excitatory and inhibitory connections. These connection indicators are multiplied by `g_ext = 0.03` and `g_inh = 0.01`, respectively, when the synaptic input is calculated.
