%% Copyright to Geetha Chandrasekaran 2023
% Author  : Geetha Chandrasekaran
% email   : geethac@utexas.edu
% Citation: doi: 10.1109/ICC51166.2024.10622169

% clear all;
% MAX_RBS = 100; %SET THIS VARIABLE TO RUN THIS FILE DIRECTLY
N_TB = 1000; % Transport block size to reflect the correct tx error prob

n_iterations = 10^4; % Increase iterations to 10^6 for higher confidence interval
n_percentiles= 15; % No. of percentiles used to track empirical CDF

% Pre-initialization of auxiliary loop variables 
needy1_RBs = zeros(1,n_iterations);
needy2_RBs = zeros(1,n_iterations);
needy3_RBs = zeros(1,n_iterations);

%% Define QoS parameters - min thruput_i, timeline tau_i
gcd_tau = 50; % The GCD of all tau_i's (tau_i: QoS measurement window user i)
tau_i     = [50 100 50 100 50 100 50 100 50 100];% tau_i's of all users
max_tau = max(tau_i); % Longest QoS measurement window
thruput_i = 1024*1000; % Target throughput
dist_rand =  [300 320 350 400 450 500 550 600 650 680];%500*ones(1,10); % fixed distance symmetric case

n_users = length(tau_i);
K_USERS = 2; % ICC paper: Support upto two deficit users before serving MaxQtl users
%tau_i = repmat(tau_i,1,n_users); % Homogeneous case tau_i's
thruput_i = repmat(thruput_i,1,n_users); % Change for variable target throughput

%% Generate arrivals, channel rates
empPercentiles = [];
rate_Percentiles = [];
tic
for i_usr=1:n_users
    % Generate random SNR for each user, and the associate data rate and Tx
    % error probabilities
    [SNR_vals(i_usr,:) rate_vals(i_usr,:) error_vals(i_usr,:)]= estimate_rates(dist_rand(i_usr), n_iterations, N_TB);
    empPercentiles = [empPercentiles; quantile(rate_vals(i_usr,:),n_percentiles)]; % Ground truth of empirical channel distribution
    [f_cdf, x_cdf]  = ecdf(rate_vals(i_usr,:));
    rate_Percentiles = [rate_Percentiles; interp1(x_cdf(2:end),f_cdf(2:end),rate_vals(i_usr,:))];
    toc
end


%% Intialize tokens needed for various policies - will be called to update after every tau_i cycle
mean_i = mean(rate_vals,2); % True mean data rates for users
Q_req_GRS = thruput_i./tau_i;
Q_req_qos = thruput_i;
Deficit = thruput_i;
serv_Deficit = thruput_i';
Power_Deficit = thruput_i';
token_queue = Q_req_GRS;
LWDF_queue = Q_req_GRS;
% 
% SLACK = max(10,floor(max_tau - ceil(n_users^2*thruput_i(1)/sum(mean(rate_vals,2)))/MAX_RBS));
% current_slack = SLACK;
% if SLACK < 0 || SLACK > max_tau
%     error('Invalid slack value')
% end

% Denominator for QoS-PF aware scheduling 
weight_den = (tau_i.*thruput_i)';

% Each user has it's own time window  tau_i  
qos_pf = zeros(1,n_iterations);

% maxRate schedule
max_rates = max(rate_vals);

avg_rates = zeros(n_users,n_iterations);
user_rates = zeros(n_users,n_iterations); % Records the throughput of scheduled user

M_lin_opt = zeros(n_users,n_iterations);
M_power = zeros(n_users,n_iterations);
M_qos = zeros(n_users,n_iterations);
M_EXP = zeros(n_users,n_iterations);
M_LWDF = zeros(n_users,n_iterations);
M_GRS = zeros(n_users,n_iterations);
M_servDeficit = zeros(n_users,n_iterations);
OGRS_def_track = zeros(n_users,n_iterations);
maxQtl_vals = zeros(1,n_iterations);


%% Oracle-aided algorithm
max_tau = max(tau_i);
for idx = 1:max_tau:n_iterations-max_tau
    % Linear programming results
    frac_mat = func_linprog(rate_vals(:,idx:idx-1+max_tau), thruput_i/MAX_RBS);
    M_lin_opt(:,idx:idx-1+max_tau) = MAX_RBS*frac_mat; 
end


%% Main routine with QoS-PF, Power-of-Two and proposed ICC paper algorithms
for idx = 1:n_iterations-max_tau
%     if idx > 100
%         [F_users(i_usr,:), x_users(i_usr,:)]  = ecdf(rate_vals(i_usr,idx-100:idx));
%     end
    OGRS_def_track(:,idx) = serv_Deficit;
    chk = repmat(rate_vals(:, idx),1,n_percentiles) - empPercentiles;
    chk = sum(sign(chk),2);
    idx_set = find(max(chk)==chk);
    mxQtl_id = idx_set(randi(length(idx_set)));    % Max Quantile id    
    maxQtl_vals(idx) = mxQtl_id;

    idx_set = find(max_rates(idx)==rate_vals(:,idx));
    mxRt_id = idx_set(randi(length(idx_set)));     % Max Rate id

    T = avg_rates(:,idx); % windowed (time) average throughput
    user_weights = rate_vals(:,idx)./ (T./weight_den+1);
    user_weights = user_weights(:).*(Deficit(:)>0);

    %% Modified PF - QoS
    % --------------------------------------------------------------------
    [~,sortIndex] = sort(user_weights,'descend');
    temp = 0;
    for s_idx = 1:n_users
        u_idx = sortIndex(s_idx) ;

        if 0 == Deficit(u_idx) % No deficit for remaining users
            break;
        end
        req_RBs = (Deficit(u_idx)/rate_vals(u_idx, idx));

        if req_RBs+temp <= MAX_RBS
            Deficit(u_idx) = 0;   % Update user service for this user
            M_qos(u_idx,idx) = (req_RBs) ;
            temp = temp + req_RBs;
        else
            M_qos(u_idx,idx) = MAX_RBS - temp ;
            Deficit(u_idx) = Deficit(u_idx) - (M_qos(u_idx,idx)*rate_vals(u_idx, idx));
            break;
        end
    end
    avg_rates(:,idx+1) = (1-1./tau_i').*avg_rates(:,idx)+1./tau_i'.*(M_qos(:,idx).*rate_vals(:, idx));

    %% Power of Two rule - for heterogeneous QoS constraints
    % --------------------------------------------------------------------
    needy_users = find(Power_Deficit);
    temp_total = 0;
    user_counter=0;
    if needy_users
        while  (sum(Power_Deficit) > 0) && (temp_total < MAX_RBS) % Assign remaining RBs to max rate
            if length(needy_users)>1
                deficit_weight = Power_Deficit; %
                deficit_weight = deficit_weight/norm(deficit_weight); % Deficit calculated for all users
                user_weights = rate_vals(needy_users,idx); % maxquantile rule
             
                idx_set = find(max(user_weights)==user_weights);
                mx_user = idx_set(randi(length(idx_set))); 
                usidx = needy_users(mx_user);
            else
                usidx = needy_users;
            end
            req_RBs = Power_Deficit(usidx)/rate_vals(usidx, idx);
            
            M_power(usidx,idx) = M_power(usidx,idx) + (min(MAX_RBS-temp_total, req_RBs)) ;
            Power_Deficit(usidx) = Power_Deficit(usidx)- (min(MAX_RBS-temp_total, req_RBs)*rate_vals(usidx,idx));

            needy_users = find(Power_Deficit);
            temp_total = sum(M_power(:,idx));
    
            %% New code bit
            user_counter = user_counter+1;
            if K_USERS < user_counter 
                break
            end
        end
    else
        % Find max rate user
        M_power(mxQtl_id,idx) = MAX_RBS ; % No slack for any user to update
    end
    % current_slack = floor((MAX_RBS*tau_i(1) - ceil(sum(Slack_Deficit)/min(rate_vals(1,:))))/MAX_RBS);
    % --------------------------------------------------------------------
    %% EXP rule
    % --------------------------------------------------------------------
    a_i = rate_vals(:,idx)./mean_i;
    avg_a_i = mean(a_i.*token_queue) ;
    user_weights = a_i.*exp((a_i.*token_queue - avg_a_i)/(1+sqrt(avg_a_i)));

    exp_total = 0;
    [~,sorted_idx] = sort(user_weights,'descend');
    for exp_idx = 1:n_users
        usr_id = sorted_idx(exp_idx);
        req_resources = token_queue(usr_id)/rate_vals(usr_id,idx);

        if req_resources + exp_total < MAX_RBS
            M_EXP(usr_id,idx) = req_resources;
            token_queue(usr_id) = 0;     
            exp_total = exp_total + req_resources;
        % Remove tokens from user's token queue 
        else
            M_EXP(usr_id,idx) = MAX_RBS - exp_total;
            token_queue(usr_id) = max(token_queue(usr_id) - M_EXP(usr_id,idx)*rate_vals(usr_id,idx), 0);
            break;
        end
        avg_a_i = mean(a_i.*token_queue) ;
        user_weights = a_i.*exp((a_i.*token_queue - avg_a_i)/(1+sqrt(avg_a_i)));
    end
    
    % Increment  the token queues every time slot
    token_queue = token_queue + Q_req_GRS;

    %% --------------------------------------------------------------------
    %% LWDF rule
    a_i = rate_vals(:,idx)./mean_i;
    user_weights = a_i.*LWDF_queue;

    % TODO: Calculate resources required
    exp_total = 0;
    [~,sorted_idx] = sort(user_weights,'descend');
    for exp_idx = 1:n_users
        usr_id = sorted_idx(exp_idx);
        req_resources = LWDF_queue(usr_id)/rate_vals(usr_id,idx);

        if req_resources + exp_total < MAX_RBS
            M_LWDF(usr_id,idx) = req_resources;
            LWDF_queue(usr_id) = 0;     
            exp_total = exp_total + req_resources;
        % Remove tokens from user's token queue 
        else
            M_LWDF(usr_id,idx) = MAX_RBS - exp_total;
            LWDF_queue(usr_id) = max(LWDF_queue(usr_id) - M_LWDF(usr_id,idx)*rate_vals(usr_id,idx), 0);
            break;
        end
        user_weights = a_i.*LWDF_queue;
    end

    % Increment  the token queues every time slot
    LWDF_queue = LWDF_queue + Q_req_GRS;

    %% Calculate GRS algo required resources
    M_GRS(:,idx) = (Q_req_GRS'./rate_vals(:,idx)) ;


    %% --------------------------------------------------------------------
    %% ICC proposed algorithm: Deficit + MaxQtl resource requirement
    temp_total = sum(M_servDeficit(:,idx));
    user_counter=0;
    needy_users = find(serv_Deficit); % Find all indices of needy users

    % IMPORTANT: If needy users, calculate no.of RBs to be allocated
    % If NO needy user, let assigned RBs be zero - max_OGRS will assign RBs
    % to user with the MaxQtl instantaneous datarate
    if needy_users
        while  (sum(serv_Deficit) > 0) && (temp_total < MAX_RBS) % Assign remaining RBs to max rate
            chk = repmat(rate_vals(needy_users, idx),1,n_percentiles) - empPercentiles(needy_users,:);
            [~, mx_user] = max(sum(sign(chk),2)); % pick maxquantile rule user idx
            usidx = needy_users(mx_user); % Update set of needy users
            req_RBs = serv_Deficit(usidx)/rate_vals(usidx, idx); % Determine RBs needed to meet deficit
            
            % Assign resources to the user selected
            M_servDeficit(usidx,idx) = M_servDeficit(usidx,idx) + (min(MAX_RBS-temp_total, req_RBs)) ;
            % Update user's deficit to reflect current RB allocation
            serv_Deficit(usidx) = serv_Deficit(usidx)- (min(MAX_RBS-temp_total, req_RBs)*rate_vals(usidx,idx));

            % Check the number of users assigned RBs
            if 1==user_counter
                needy2_RBs(idx) = min(MAX_RBS-temp_total, req_RBs);
            elseif 0==user_counter
                needy1_RBs(idx) = min(MAX_RBS-temp_total, req_RBs);
            else
                needy3_RBs(idx) = min(MAX_RBS-temp_total, req_RBs);
            end
            
            % Update the set of needy users
            needy_users = find(serv_Deficit);
            temp_total = sum(M_servDeficit(:,idx));
    
            % Update number of needy users supported in this iteration
            user_counter = user_counter+1;
            if K_USERS < user_counter 
                break
            end
        end
    end

    % If at the end of tau_i additional Deficit
    if 0 == mod(idx,gcd_tau)
        tau_cycle = mod(idx,tau_i);
        for cycle_idx = 1:n_users
            % Check the token buckets that need to be re-initialized based
            % on tau_i cycles
            if 0 == tau_cycle(cycle_idx)
                Power_Deficit =  Q_req_qos(cycle_idx);
                Deficit(cycle_idx) = Q_req_qos(cycle_idx);
                serv_Deficit(cycle_idx) = Q_req_qos(cycle_idx);
            end
        end
    end
end

    
thruput_GRS = M_GRS.*rate_vals;
thruput_OGRS = M_servDeficit.*rate_vals;
thruput_power = M_power.*rate_vals;
thruput_EXP = M_EXP.*rate_vals;
thruput_LWDF = M_LWDF.*rate_vals;
thruput_QoS = M_qos.*rate_vals;
thruput_lin_opt = M_lin_opt.*rate_vals;

% Find max quantile users for the last few iterations
for idx = n_iterations-max_tau+1:n_iterations
    OGRS_def_track(:,idx) = serv_Deficit;
    chk = repmat(rate_vals(:, idx),1,n_percentiles) - empPercentiles;
    chk = sum(sign(chk),2);
    idx_set = find(max(chk)==chk);
    mxQtl_id = idx_set(randi(length(idx_set)));    % Max Quantile id    
    maxQtl_vals(idx) = mxQtl_id;
end 

figure
ecdf(sum(M_GRS)); hold all;
ecdf(sum(M_servDeficit));
ecdf(sum(M_qos));
title('CDF of resource requirement')

figure
plot(cumsum(sum(M_GRS)),'DisplayName','GRS')
hold all;
plot(cumsum(sum(M_servDeficit)),'DisplayName','OGRS')
plot(cumsum(sum(M_qos)),'DisplayName','QoS-PF');
plot(cumsum(sum(M_lin_opt)),'DisplayName','LinProg');  hold off;
legend

max_OGRS = (MAX_RBS - sum(M_servDeficit)).*max(rate_vals);
max_QoS = (MAX_RBS - sum(M_qos)).*max(rate_vals);
max_EXP = (MAX_RBS - sum(M_EXP)).*max(rate_vals);
max_LWDF = (MAX_RBS - sum(M_LWDF)).*max(rate_vals);
max_linOpt = (MAX_RBS - sum(M_lin_opt)).*max(rate_vals);
max_power = (MAX_RBS - sum(M_power)).*max(rate_vals);

max_tpt_OGRS = zeros(size(M_servDeficit));
max_tpt_QoS = zeros(size(M_servDeficit));
max_tpt_opt = zeros(size(M_servDeficit));
max_tpt_power = zeros(size(M_servDeficit));
max_tpt_EXP = zeros(size(M_servDeficit));
max_tpt_LWDF = zeros(size(M_servDeficit));
maxrate_usrs = zeros(size(M_servDeficit));

% Set the user rate correctly
% for i = 1:n_iterations
%     idx_set = find(max_rates(i)==rate_vals(:,i));
%     mx_id = idx_set(randi(length(idx_set)));
%     max_tpt_QoS(mx_id,i) = max_QoS(i);
%     max_tpt_opt(mx_id,i) = max_linOpt(i);
%     max_tpt_OGRS(mx_id,i) = max_OGRS(i);
%     max_tpt_EXP(mx_id,i) = max_EXP(i);
%     max_tpt_LWDF(mx_id,i) = max_LWDF(i);
%     max_tpt_power(mx_id,i) = max_power(i);
%     maxrate_usrs(mx_id,i) = MAX_RBS;
% end

% Resolve user ties correctly
for i = 1:n_iterations
    max_tpt_QoS(maxQtl_vals(i),i) = max_QoS(i);
    max_tpt_opt(maxQtl_vals(i),i) = max_linOpt(i);
    max_tpt_OGRS(maxQtl_vals(i),i) = max_OGRS(i); %maxQtl
    max_tpt_EXP(maxQtl_vals(i),i) = max_EXP(i); %maxQtl
    max_tpt_LWDF(maxQtl_vals(i),i) = max_LWDF(i); %maxQtl
    max_tpt_power(maxQtl_vals(i),i) = max_power(i); %maxQtl
    maxrate_usrs(maxQtl_vals(i),i) = MAX_RBS; %maxQtl
end

% Calculate performance metrics
for i = 1:n_users
    perform_GRS{i,:} = sum(reshape(thruput_GRS(i,:), ...
        [tau_i(i), n_iterations/tau_i(i)]));
    perform_OGRS{i,:} = sum(reshape(thruput_OGRS(i,:)+max_tpt_OGRS(i,:), ...
        [tau_i(i), n_iterations/tau_i(i)]));
    perform_QoS{i,:} = sum(reshape(thruput_QoS(i,:)+max_tpt_QoS(i,:), ...
        [tau_i(i), n_iterations/tau_i(i)]));
    perform_power{i,:} = sum(reshape(thruput_power(i,:)+max_tpt_power(i,:), ...
        [tau_i(i), n_iterations/tau_i(i)]));
    perform_maxrate{i,:} = sum(reshape(maxrate_usrs(i,:), ...
        [tau_i(i), n_iterations/tau_i(i)]));
    perform_EXP{i,:} = sum(reshape(thruput_EXP(i,:)+ max_tpt_EXP(i,:), ...
        [tau_i(i), n_iterations/tau_i(i)]));
    perform_LWDF{i,:} = sum(reshape(thruput_LWDF(i,:)+ max_tpt_LWDF(i,:), ...
        [tau_i(i), n_iterations/tau_i(i)]));
end

% Calculate the mean rates allocated by each algorithm
mean_OGRS = mean(thruput_OGRS+max_tpt_OGRS,2); %Best heterogeneous algo
mean_QoS = mean(thruput_QoS+max_tpt_QoS,2);
mean_lin_opt = mean(thruput_lin_opt+max_tpt_opt,2);
mean_power = mean(thruput_power+max_tpt_power,2);
mean_EXP = mean(thruput_EXP+max_tpt_EXP,2);
mean_LWDF = mean(thruput_LWDF+max_tpt_LWDF,2);
target_rate = mean(thruput_GRS,'all');
sum(mean_OGRS)
close all;

% Create and update mean rate CDF plots across users
[f,x] = ecdf(mean_EXP); plot(x,f,'LineStyle','-','Marker','o');  hold all;
[f,x] = ecdf(mean_LWDF); plot(x,f,'LineStyle','-','Marker','o');  hold all;
[f,x] = ecdf(mean_QoS); plot(x,f,'LineStyle','-','Marker','^'); 
[f,x] = ecdf(mean_OGRS); plot(x,f,'LineStyle','-','Marker','+');
ecdf(mean_power); 
[f,x] = ecdf(mean(maxrate_usrs.*rate_vals,2)); plot(x,f,'LineStyle',':'); 
[f,x] = ecdf(mean_lin_opt);  plot(x,f,'LineStyle','--'); 
legend('LWDF','QoS-PF','Deficit','Power_of_2','MaxQtl','Oracle','Location', 'Best')
xlabel('Mean user Rate')
ylabel('Empirical CDF')
box off;
xlim([0 6*10^4])

% Calculate the number of Tx errors across users
error_EXP = [];
error_power = [];
error_QoS_PF = [];
error_OGRS = [];
error_LWDF = [];
for i = 1:n_users
    error_LWDF = [error_LWDF sum(perform_LWDF{i,:}<thruput_i(i)-0.01)/n_iterations*tau_i(i)];
    error_EXP = [error_EXP sum(perform_EXP{i,:}<thruput_i(i)-0.01)/n_iterations*tau_i(i)];
    error_power = [error_power  sum(perform_power{i,:}<thruput_i(i)-0.01)/n_iterations*tau_i(i)];
    error_QoS_PF = [error_QoS_PF   sum(perform_QoS{i,:}<thruput_i(i)-0.01)/n_iterations*tau_i(i)];
    error_OGRS = [error_OGRS sum(perform_OGRS{i,:}<thruput_i(i)-0.01)/n_iterations*tau_i(i)];
end
%save(strcat('error10_7_tau',num2str(max_tau),num2str(MAX_RBS),'_maxRBs'),'error_*');
%save(strcat('Trial_10_6_USERS_',num2str(K_USERS),'.mat'),'error*'); %strcat('TrueHet_err10_7_',num2str(MAX_RBS),'.mat'),'error*');
error_OGRS
toc

% figure (DEBUG ONLY)
% need_3  = find(needy3_RBs~=0);
% Y = [needy1_RBs' needy2_RBs' needy3_RBs' (MAX_RBS-needy1_RBs-needy2_RBs-needy3_RBs)'];
% bar(M_OGRS(:,1:100)','stacked')

% Function that provides the Oracle-aided algorithm to satisfy QoS of all
% users
function frac_mat = func_linprog(ch_mat, thruput_i)
    % Optimization routine for given time window and min rate requirement
    n_users = size(ch_mat,1);
    time_wdw = size(ch_mat,2);
    ch_vec = ch_mat(:);
    func = ch_vec';
    
    A_rep = [ones(1,n_users) repmat(zeros(1,n_users),1,(time_wdw-1))];
    A = A_rep;
    for j=1:time_wdw-1
        A = [A; circshift(A_rep,[0 n_users*j])];
    end
    
    B_rep = repmat([1 zeros(1,n_users-1)],1,time_wdw);
    B = B_rep.*ch_vec';
    for j=1:n_users-1
        B_rep = circshift(B_rep,1);
        B = [B; B_rep.*ch_vec'];
    end
    A_con = [A; -B];
    b_con = [ones(time_wdw,1); -thruput_i(:)]; % zero lower bound
    lb = zeros(1,n_users*time_wdw); 
    ub = lb+1; %one upper bound
    options = optimoptions('linprog','Display','none');
    x = linprog(-func,A_con,b_con,[],[],lb,ub,options);
    try
        frac_mat = reshape(x,size(ch_mat));
    catch
        disp('Optimal solution not found')
        frac_mat = zeros(size(ch_mat));
    end
    
end
