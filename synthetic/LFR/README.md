# How to generate the LFR Benchmark graphs

1. Clone the repository [LFRbenchmarks](https://github.com/andrealancichinetti/LFRbenchmarks/tree/main).
1. Enter in the directory `unweighted_directed`.
1. Copy the script `generate_LFR.sh` inside the folder and make it executable:
```bash
chmod +x generate_LFR.sh
```
1. Run the script:
```bash
./generate_LFR.sh
```
1. Move all the generated graphs from the folder `graphs` to the folder `synthetic/LFR/raw_graphs` of the spectral-clustering repository.
1. Run the MATLAB script `get_LFR_graphs.m`.
1. All the graphs compatible with our procedure will be available inside the folder `synthetic/LFR/generated`.
