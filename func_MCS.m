%% Copyright information
% Author : Geetha Chandrasekaran
% email  : geethac@utexas.edu
% Website: https://scholar.google.com/citations?user=kOI1ZGkAAAAJ
% Last revision: Oct 13, 2023.
% Add citation: doi: 10.1109/ICC51166.2024.10622169
% G. Chandrasekaran and G. d. Veciana, "Opportunistic Scheduling 
% for Users with Heterogeneous Minimum Rate QoS Requirements," 
% ICC 2024 - IEEE International Conference on Communications, Denver, CO, USA, 2024, pp. 1-6. 

%% Coding rate and MCS mapping based on SNR
% -----------------------------------------
% Input parameters:
% -------------------------------------------------------------------------
% 1. SNR in linear scale
% 2. b3GPP = 1 (3GPP pathloss), otherwise linear interpolated
% 3. NTB = transport block size
% -------------------------------------------------------------------------
% Output parameters:
% quantized SNR, Quantized rate, Polyanski error probability
% -----------------------------------------
% MCS index based on the follwing tables
%
%---------------------------------------
%Index Mod    Coding(bits/sym) SNR (dB)
%---------------------------------------
%  0  QPSK    -    0        - \\
%  1  QPSK    78   0.1523   -6.934 \\
%  2  QPSK    120  0.2344   -5.147 \\
%  3  QPSK    193  0.3770   -3.18 \\
%  4  QPSK    308  0.6016   -1.254\\
%  5  QPSK    409  0.8770    0.761 \\
%  6  QPSK    602  1.1758    2.697 \\
%  7  16-QAM  378  1.4766    4.697  \\
%  8  16-QAM  490  1.9141    6.528\\
%  9  16-QAM  616  2.4063    8.576 \\ 
% 10  64-QAM  466  2.7305    10.37 \\
% 11  64-QAM  567  3.3223    12.3 \\
% 12  64-QAM  666  3.9023    14.18 \\
% 13  64-QAM  772  4.5234    15.89\\
% 14  64-QAM  872  5.1152    17.82\\
% 15  64-QAM  948  5.5547    19.83\\

function [snr, rate, p_err] = func_MCS(SNR, b3gpp, N_TB) % 
    bits_per_symbol= [ 0.1523/8
    0.1523/4
    0.1523/2
    0.1523
    0.2344
    0.3770
    0.6016
    0.8770
    1.1758
    1.4766
    1.9141
    2.4063
    2.7305
    3.3223
    3.9023
    4.5234
    5.1152
    5.5547];

    mcs_vals = 1:18;
    if b3gpp %3GPP pathloss model
        snr_vals = 10.^([-20 -15 -10 -6.934 -5.147 -3.18 -1.254 0.761 ...
            2.697 4.697 6.528 8.576 10.37 ...
            12.3 14.18 15.89 17.82 19.83 ]/10);
        rate_vals = [18 36 2*36 2*78 2*120 2*193 2*308 2*409 ...
            2*602 4*378 4*490 4*616 6*466 ...
            6*567 6*666 6*772 6*872 6*948];
    else
        snr_vals = linspace(10^(-6.934 /10), 100, length(mcs_vals));
        rate_vals = log2(1+snr_vals);
    end
    
    snr_index = find(SNR-snr_vals>=0);
    b_sym =  bits_per_symbol(snr_index(end));

    % Find the probability of error as a function of SNR (Polyanski)
    q_val = (log2(1+SNR) - b_sym)/(log2(exp(1))*sqrt((1-1/(1+SNR)^2)/N_TB));
    p_err = qfunc(q_val);
    
    if snr_index  
        mcs = mcs_vals(snr_index(end));
        snr = snr_vals(snr_index(end));
        rate = rate_vals(mcs);
        
    else
        snr = 0;
        rate = 0;
    end

end
