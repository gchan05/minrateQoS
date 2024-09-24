%% Copyright information
% Author : Geetha Chandrasekaran
% email  : geethac@utexas.edu
% Website: https://scholar.google.com/citations?user=kOI1ZGkAAAAJ
% Last revision: Oct 13, 2023.
% Add citation: doi: 10.1109/ICC51166.2024.10622169
% G. Chandrasekaran and G. d. Veciana, "Opportunistic Scheduling 
% for Users with Heterogeneous Minimum Rate QoS Requirements," 
% ICC 2024 - IEEE International Conference on Communications, Denver, CO, USA, 2024, pp. 1-6. 

%% Function to generate 3GPP pathloss or standard exponential pathloss
% d_2D = 2D distance between the eNB and the UE
% d_3D = 3D distance between the eNB and the UE
% fc = carrier frequency
% hBS(hUE) = height of BS(user equipment)
% -----------------------------------------
% Output parameters:
% -----------------------------------------
% pathloss in dB
%%
function [pathloss, isLOS]= func_PathLoss(d_2D,d_3D,fc,hBS,hUE,functionFlag)

alpha = 4; % Pathloss exponent

if strcmp(functionFlag,'3GPP')
    N_BS=size(d_2D,1);
    N_UE=size(d_2D,2);
    
    % Distance threshold
    d_BP = 4*(hBS-1)*(hUE-1)*fc*10^9/(3*10^8); %use heights in m, fc in Hz
    probLOS = 18./d_2D+exp(-d_2D/63).*(1-18./d_2D);
    
    isLOS = rand(N_BS,N_UE)<probLOS;
    closeThan_d_BP = d_2D<=d_BP;
    Further_d_BP = ~closeThan_d_BP;
    
    % LoS if close to the BS
    isLOS(closeThan_d_BP)= 1;
    
    PL1 = 32.4+21*log10(d_3D)+20*log10(fc);
    PL2 = 32.4+40*log10(d_3D)+20*log10(fc)-9.5*log10(d_BP^2+(hBS-hUE)^2);
    
    % Calculate LOS pathloss
    PL_LOS = PL1; %Set to PL1 for initialization
    PL_LOS(Further_d_BP) = PL2(Further_d_BP);
    
    % Calculate NLOS pathloss
    PL3 = 35.3*log10(d_3D)+22.4+21.3*log10(fc)-0.3*(hUE-1.5);
    PL_NLOS = max(PL_LOS,PL3);
    
    % Path loss in dB
    pathloss = PL_LOS.*isLOS+PL_NLOS.*(1-isLOS);
    
elseif strcmp(functionFlag,'standard')
    % Standard path loss in dB
    pathloss= alpha*10*(log10(d_3D)+log10(fc)+log10(4*pi/3e8));
else
    error('Undefined Pathloss Type');
    
end


end