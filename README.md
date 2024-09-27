# Rate QoS constrained Scheduling
This repository contains all necessary files required to generate the plots provided in G. Chandrasekaran and G. d. Veciana, "Opportunistic Scheduling for Users with Heterogeneous Minimum Rate QoS Requirements," ICC 2024 - IEEE International Conference on Communications, Denver, CO, USA, 2024, pp. 1-6, doi: 10.1109/ICC51166.2024.10622169.

Matlab 2023b was used to generate the attached plot, for details on the simulation parameters and scheduler benchmarks please refer to the above ICC conference paper.

# Run instructions
1. Download the entire folder
2. Run the file eMBB_main.m to generate plots
   
I recommend setting the variable n_iterations = 10^4 ('Hetro_eMBB.m', Line 10) to a much higher value (like 10^6) to obtain results with a high confidence interval.
