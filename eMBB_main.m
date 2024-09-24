%% Copyright information
% Author : Geetha Chandrasekaran
% email  : geethac@utexas.edu
% Website: https://scholar.google.com/citations?user=kOI1ZGkAAAAJ
% Last revision: Oct 13, 2023.
% Add citation: doi: 10.1109/ICC51166.2024.10622169
% G. Chandrasekaran and G. d. Veciana, "Opportunistic Scheduling 
% for Users with Heterogeneous Minimum Rate QoS Requirements," 
% ICC 2024 - IEEE International Conference on Communications, Denver, CO, USA, 2024, pp. 1-6. 

MAX_RB_vals = 60;%80:20:300;

for i = 1:length(MAX_RB_vals)
    MAX_RBS = MAX_RB_vals(i);
    Hetro_eMBB
    %Hetro_variations
end