# Heuristics-and-metaheuristics-for-Acoustic-Topology-Optimisation

This repository contains the source code for the research work titled "Comparison of heuristics and metaheuristics for topology optimisation in acoustic porous materials" published in the Journal of Acoustical Society of America (JASA).

### Authors: 
- Vivek T Ramamoorthy
- Ender Özcan
- Andrew J. Parkes
- Abhilash Sreekumar
- Luc Jaouen
- François-Xavier Bécot

### DOI

[https://doi.org/10.1121/10.0006784](https://doi.org/10.1121/10.0006784)



## Documentation
The matlab code main files are organised as follows:

### Main file:

`benchmarktest.m` is the starting-point script to run the optimisation trials. 
Running this after uncommenting the `run_<algorithm_name>` should run and store results while also plotting them.


### Other files: 
`benchmarks.m` is a function that returns the problem instance information from a given problem instance number.

### FE Solvers:
`AFSO_input_from_benchmarks.m` is a script that uses the instance information from `benchmarks.m` and sets up the necessary variables to call the solver


`AFSO_solver_BH_comb.m` is the Biot-Hemlholtz combinatorial solver function that computes and returns the frequency-averaged absorption values given an input shape. This only accepts 0 or 1 inputs for the shape


`AFSO_solver_BH_cont.m` is the corresponding continous solver function that accepts shapes with between 0 and 1 design variables. Also, if three outputs are requested, then the third output is the gradient, the second output is the absorption vector at all frequencies

### Decoders or Solver wrappers:

`griddecoder_comb.m` and `griddecoder_cont.m` are functions that take the shape in terms of a binary array as input and nothing else, and return the frequency-averaged sound absorption. The latter also returns the absorption array at all frequencies as the second output, and the sensitivities of absorption as the third output. Note the the derivative computation takes significant time and the code automatically skips the computation if only one output is requested.

### Algorithm runners
All files with the prefix `run_` are algorithm running scripts that run trials and save them in a heirarchical folder structure while also doing post processing and plotting the results. It is advised to call these scripts from the `benchmarktest.m` file rather than standalone.


## Acknowledgement
This work was supported by the European Union's Horizon 2020 Marie Sklodowska Curie Actions grant 765472 ITN No2Noise ([no2noise.eu](https://no2noise.eu))

## License
GPL-3.0 License



