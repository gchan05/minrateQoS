%% Copyright information
% Author : Geetha Chandrasekaran
% email  : geethac@utexas.edu
% Website: https://scholar.google.com/citations?user=kOI1ZGkAAAAJ
% Last revision: Oct 13, 2023.
% Add citation: doi: 10.1109/ICC51166.2024.10622169
% G. Chandrasekaran and G. d. Veciana, "Opportunistic Scheduling 
% for Users with Heterogeneous Minimum Rate QoS Requirements," 
% ICC 2024 - IEEE International Conference on Communications, Denver, CO, USA, 2024, pp. 1-6. 

function [SNR_vals, rate_vals, error_vals] = estimate_rates(d_2D, nIter, PKT_SIZE)

    txPwr = 2; %BS transmit power 2W
    BS_HT = 10; %5m
    UE_HT = 1.6; %1.6m
    NOISE_PWR = 10^(-10.4)*10^(-3); % -104 dBm white noise power
    fc = 3.5; % 3.5GHz;
    BS_loc = [0 0 BS_HT];
    t = 1:nIter;%;10^2/2; % time index
    ch_vals = ones(size(t));
    rate_vals = ones(size(t));
    error_vals = ones(size(t));
    
    %% UE location
    theta = rand*2*pi; %random angle between (0,2*pi)
    UE_loc = [ d_2D*cos(theta)  d_2D*sin(theta) UE_HT];
    d_3D = sqrt(sum((BS_loc-UE_loc).^2));

    % generate 3GPP random pathloss and small scale fading
    pl_vals = func_PathLoss(repelem(d_2D,1,10*nIter), repelem(d_3D, 1,10*nIter), fc, BS_HT, UE_HT, '3GPP');
    fading_vals = exprnd(1,1,10*nIter) ; %|h|^2 ~ Exp(1) for Rayleigh fading
    SNR_vals = (txPwr*10.^(-pl_vals/10).*fading_vals/NOISE_PWR);
    
    SNR_MIN = 10^(-6.93/10); % linear scale
    SNR_MAX  = 10^(19.83/10); % linear scale
    
    % bounding the SNR values between SNR_MIN and SNR_MAX
    SNR_vals = SNR_vals(SNR_vals>SNR_MIN);
    SNR_vals = SNR_vals(1:nIter);
    SNR_vals(SNR_vals>SNR_MAX) = SNR_MAX;
    %SNR_normalized = rescale(SNR_vals,'InputMin',10^(-0.697),'InputMax',10^(2));
    
    for i=1:length(SNR_vals)
        %ch_vals(i) = SNR_vals(i);
        %rate_vals(i) = log2(1+SNR_vals(i));
        SNR = SNR_vals(i);
        [ch_vals(i), rate_vals(i), error_vals(i)] = func_MCS(SNR,1, PKT_SIZE); %TODO: change for 3GPP
    end                                                                                                                                                     
end