%% Copyright information
% Author : Geetha Chandrasekaran
% email  : geethac@utexas.edu
% Website: https://scholar.google.com/citations?user=kOI1ZGkAAAAJ
% Last revision: Oct 13, 2023.
% Add citation: doi: 10.1109/ICC51166.2024.10622169
% G. Chandrasekaran and G. d. Veciana, "Opportunistic Scheduling 
% for Users with Heterogeneous Minimum Rate QoS Requirements," 
% ICC 2024 - IEEE International Conference on Communications, Denver, CO, USA, 2024, pp. 1-6. 

%% TODO: Can change arrival pattern to real time video traffic
% clear all;
n_iterations = 10^3;
n_percentiles= 15;

MAX_RBS = 300;

gcd_tau = 50;
tau_i     = 50;
thruput_i = 1024*1000;
dist_rand = 500*ones(1,5); % fixed distance symmetric case

n_users = 10;
tau_i = repmat(tau_i,1,n_users);
thruput_i = repmat(thruput_i,1,n_users);
dist_rand = repmat(dist_rand,1,n_users);

%% Generate arrivals, channel rates
empPercentiles = [];
tic
for i_usr=1:n_users
    %[rate_vals(i_usr,:)]= estimate_rates(dist_rand(i_usr), n_iterations);
    empPercentiles = [empPercentiles; quantile(rate_vals(i_usr,:),n_percentiles)];
    toc
end

%% Define QoS parameters - min thruput_i, timeline tau_i
Q_req_GRS = thruput_i./tau_i;
a_i = thruput_i./tau_i;
Q_req_qos = thruput_i;
Deficit = thruput_i;
OGRS_Deficit = thruput_i';
Slack_Deficit = thruput_i';
token_queue = Q_req_GRS;
schedule_wdw = 10;

SLACK = floor((MAX_RBS*tau_i(1) - n_users*ceil(thruput_i(1)/mean(rate_vals(1,:))))/MAX_RBS);
current_slack = SLACK;
if SLACK < 0
    error('Invalid slack value')
end

% Denominator for QoS aware scheduling 
weight_den = (tau_i.*thruput_i)';

% Each user has it's own time window  tau_i  
qos_pf = zeros(1,n_iterations);

% maxRate schedule
max_rates = max(rate_vals);

avg_rates = zeros(n_users,n_iterations);
user_rates = zeros(n_users,n_iterations); % Records the throughput of scheduled user

M_heuristic = zeros(n_users,n_iterations);
M_lin_opt = zeros(n_users,n_iterations);
M_slack = zeros(n_users,n_iterations);
M_qos = zeros(n_users,n_iterations);
M_GRS = zeros(n_users,n_iterations);
M_OGRS = zeros(n_users,n_iterations);
OGRS_def_track = zeros(n_users,n_iterations);

for idx = 1:n_iterations-tau_i(1)
    OGRS_def_track(:,idx) = OGRS_Deficit;
    chk = repmat(rate_vals(:, idx),1,n_percentiles) - empPercentiles;
    [~, mxQtl_user] = max(sum(sign(chk),2)); % maxquantile rule

    T = avg_rates(:,idx); % windowed (time) average throughput
    user_weights = rate_vals(:,idx)./ (T./weight_den+1);
    user_weights = user_weights(:).*(Deficit(:)>0);

    %% Modified PF - QoS
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

    %% SLACK rule - TODO:? change for heterogeneous case
    if idx == 620
        disp('Point of trouble');
    end
    if mod(idx,tau_i(1)) < SLACK
        idx_set = find(max_rates(idx)==rate_vals(:,idx));
        mx_id = idx_set(randi(length(idx_set)));
        M_slack(mx_id,idx) = MAX_RBS ;
        % update slack if required
        if Slack_Deficit(mx_id)
            Slack_Deficit(mx_id) = max(Slack_Deficit(mx_id)-MAX_RBS*rate_vals(mx_id,idx),0);
        end
    else
        needy_users = find(Slack_Deficit);
        if needy_users
            temp_total = 0;
            % Find maxRate user and allocate RBs - then fill up remaining RBs
            [~,sorted_idx] = sort(rate_vals(needy_users,idx)+Slack_Deficit(needy_users)./thruput_i(needy_users),'descend'); %% TODO: test this piece of code
            needy_sorted = needy_users(sorted_idx);
            for temp_idx = 1:length(needy_sorted)
                mx_id = needy_sorted(temp_idx);
                req_resources = Slack_Deficit(mx_id)/rate_vals(mx_id,idx);
                if req_resources < temp_total
                    % with next best user and so on
                    M_slack(mx_id,idx) = req_resources ;
                    Slack_Deficit(mx_id) = 0; % No pendng Deficit 
        
                    temp_total = temp_total+ req_resources;
                else
                    M_slack(mx_id,idx) = MAX_RBS-temp_total ;
                    Slack_Deficit(mx_id) = max(Slack_Deficit(mx_id)-MAX_RBS*rate_vals(mx_id,idx),0);
                    break;
                end
            end %end-for temp_idx

        else
            % Find max rate user
            idx_set = find(max_rates(idx)==rate_vals(:,idx));
            mx_id = idx_set(randi(length(idx_set)));
            M_slack(mx_id,idx) = MAX_RBS ; % No slack for any user to update
        end
    end
    % current_slack = floor((MAX_RBS*tau_i(1) - ceil(sum(Slack_Deficit)/min(rate_vals(1,:))))/MAX_RBS);
    %% --------------------------------------------------------------------
    %% EXP rule
    avg_a_i = a_i.*token_queue;
    user_weights = a_i.*exp((a_i.*token_queue - avg_a_i)/(1+sqrt(avg_a_i)));

    % TODO: Calculate resources required

    % Remove tokens from user's token queue

    %% Calculate GRS algo required resources
    M_GRS(:,idx) = (Q_req_GRS'./rate_vals(:,idx)) ;


    %% --------------------------------------------------------------------
    %% Deficit + Maxrate resource requirement
    temp_total = sum(M_OGRS(:,idx));
    user_counter=0;
    while  (sum(OGRS_Deficit) > 0) && (temp_total < MAX_RBS) % Assign remaining RBs to max rate
        needy_users = find(OGRS_Deficit); % Find all indices of needy users
        chk = repmat(rate_vals(needy_users, idx),1,n_percentiles) - empPercentiles(needy_users,:);
        [~, mx_user] = max(sum(sign(chk),2)); % maxquantile rule
%        [~, mx_user] = max(rate_vals(needy_users, idx)); % maxrate rule
        usidx = needy_users(mx_user);
        req_RBs = OGRS_Deficit(usidx)/rate_vals(usidx, idx);
        
        M_OGRS(usidx,idx) = M_OGRS(usidx,idx) + (min(MAX_RBS-temp_total, req_RBs)) ;
        OGRS_Deficit(usidx) = OGRS_Deficit(usidx)- (min(MAX_RBS-temp_total, req_RBs)*rate_vals(usidx,idx));
        needy_users = find(OGRS_Deficit);
        temp_total = sum(M_OGRS(:,idx));

        %% New code bit
        user_counter = user_counter+1;
        if 2 < user_counter 
            break
        end
    end
%     [M_est_sorted, sorted_Index] = sort(M_est_OGRS);
%     last_idx = find(MAX_RBS - cumsum(M_est_sorted)<0,1);
%     if last_idx
%         for ogrs_idx = 1:last_idx-1
%             M_OGRS(sorted_Index(ogrs_idx),idx) = M_est_OGRS(sorted_Index(ogrs_idx)); %set for each user
%         end
%     else
%         M_OGRS(:,idx) = M_est_OGRS;
%     end
%    OGRS_Deficit = OGRS_Deficit - M_OGRS(:,idx)'.*rate_vals(:,idx)';

    % Increment  the token queues every time slot
    token_queue = token_queue + Q_req_GRS;

    % If at the end of tau_i additional Deficit
    if 0 == mod(idx,gcd_tau)
        tau_cycle = mod(idx,tau_i);
        for cycle_idx = 1:n_users % Reset slack values 
            if 0 == tau_cycle(cycle_idx)
                current_slack = SLACK; %%TODO: change for heterogeneous case
                Deficit(cycle_idx) =  Q_req_qos(cycle_idx);
                OGRS_Deficit(cycle_idx) = Q_req_qos(cycle_idx);
                Slack_Deficit(cycle_idx) =  Q_req_qos(cycle_idx);
            end
        end
    end
end

for idx = 1:tau_i(1):n_iterations-tau_i(1)
    % Linear programming results
    frac_mat = func_linprog(rate_vals(:,idx:idx-1+tau_i(1)), thruput_i/MAX_RBS);
    M_lin_opt(:,idx:idx-1+tau_i(1)) = MAX_RBS*frac_mat; 

    %for each user find the best slot in the next tau_(1) slots
    Opt_deficit = Q_req_qos;
    for usr_i = 1:n_users
        [B, I] = sort(rate_vals(usr_i, idx:idx-1+tau_i(1)),'descend');
        max_val = B(1); max_idx = I(1);
        total_current = MAX_RBS - sum(M_heuristic(:,idx+max_idx));
        req_RBS = Opt_deficit(usr_i)/max_val;
        if req_RBS < total_current
            M_heuristic(usr_i,idx+max_idx) = req_RBS; %Done with the requirement
            Opt_deficit(usr_i) = 0;
        else
            % Find the next best time slot available
            M_heuristic(usr_i,idx+max_idx) = (MAX_RBS- total_current)/max_val;
            Opt_deficit(usr_i) = Opt_deficit(usr_i) - max_val*M_heuristic(usr_i,idx+max_idx);
            for j = 2:length(I)
                next_slot = I(j);
                next_rate = B(j);
                total_current = sum(M_heuristic(:,idx+next_slot));
                if  total_current < MAX_RBS
                    M_heuristic(usr_i,idx+next_slot) = min(MAX_RBS- total_current, Opt_deficit(usr_i))/next_rate;
                    Opt_deficit(usr_i) = Opt_deficit(usr_i) - next_rate*M_heuristic(usr_i,idx+next_slot);
                end
                if 0 == Opt_deficit(usr_i)
                    break;
                end
            end
        end
    end
end

thruput_GRS = M_GRS.*rate_vals;
thruput_OGRS = M_OGRS.*rate_vals;
thruput_Slack = M_slack.*rate_vals;
thruput_QoS = M_qos.*rate_vals;
thruput_heu = M_heuristic.*rate_vals;
thruput_lin_opt = M_lin_opt.*rate_vals;


figure
ecdf(sum(M_GRS)); hold all;
ecdf(sum(M_OGRS));
ecdf(sum(M_qos));
title('CDF of resource requirement')

% figure
% title('Aggregate Thruput')
% ecdf(sum(thruput_GRS))
% hold on;
% ecdf(sum(thruput_OGRS));
% ecdf(sum(thruput_QoS));hold off;
% 
figure
plot(cumsum(sum(M_GRS)),'DisplayName','GRS')
hold all;
plot(cumsum(sum(M_OGRS)),'DisplayName','OGRS')
plot(cumsum(sum(M_qos)),'DisplayName','QoS');
plot(cumsum(sum(M_heuristic)),'DisplayName','Heuristic');
plot(cumsum(sum(M_lin_opt)),'DisplayName','LinProg');  hold off;
legend

%Q_qos(:,idx) = Q_qos(:,idx-1) - D_qos(:,idx-1) + A(:,idx-1) ;
% f = @(X) sum(X, 'all');
% g = @(X, H, rate, X_max) [(rate - sum(X.*H,2)) ; (sum(X) - X_max )'];
% gcon = @(X) g(X, CC, thruput_i(:), MAX_RBS);

max_OGRS = (MAX_RBS - sum(M_OGRS)).*max(rate_vals);
max_QoS = (MAX_RBS - sum(M_qos)).*max(rate_vals);
max_heuristic = (MAX_RBS - sum(M_heuristic)).*max(rate_vals);
max_slack = (MAX_RBS - sum(M_slack)).*max(rate_vals);

max_tpt_OGRS = zeros(size(M_OGRS));
max_tpt_QoS = zeros(size(M_OGRS));
max_tpt_opt = zeros(size(M_OGRS));
max_tpt_slack = zeros(size(M_OGRS));

% Set the user rate correctly
for i = 1:n_iterations
    idx_set = find(max_rates(i)==rate_vals(:,i));
    mx_id = idx_set(randi(length(idx_set)));
    max_tpt_OGRS(mx_id,i) = max_OGRS(i);
    max_tpt_QoS(mx_id,i) = max_QoS(i);
    max_tpt_opt(mx_id,i) = max_heuristic(i);
    max_tpt_slack(mx_id,i) = max_slack(i);
    maxrate_usrs(mx_id,i) = MAX_RBS;
end

% Calculate performance metrics
for i = 1:n_users

    perform_GRS(i,:) = sum(reshape(thruput_GRS(i,:), ...
        [tau_i(i), n_iterations/tau_i(i)]));
    perform_OGRS(i,:) = sum(reshape(thruput_OGRS(i,:)+max_tpt_OGRS(i,:), ...
        [tau_i(i), n_iterations/tau_i(i)]));
    perform_QoS(i,:) = sum(reshape(thruput_QoS(i,:)+max_tpt_QoS(i,:), ...
        [tau_i(i), n_iterations/tau_i(i)]));
    perform_slack(i,:) = sum(reshape(thruput_Slack(i,:)+max_tpt_slack(i,:), ...
        [tau_i(i), n_iterations/tau_i(i)]));
    perform_maxrate(i,:) = sum(reshape(maxrate_usrs(i,:), ...
        [tau_i(i), n_iterations/tau_i(i)]));

    figure(i)
    plot(perform_GRS(i,:),':','DisplayName','GRS'); hold on;
    plot(perform_OGRS(i,:),'DisplayName','OGRS'); 
    plot(perform_slack(i,:),'DisplayName','Slack'); 
    plot(perform_QoS(i,:),'--','DisplayName','QoS');
    plot(perform_maxrate(i,:),'--','DisplayName','MaxRate');
    legend;
    % hold off;
end

mean_OGRS = mean(thruput_OGRS+max_tpt_OGRS,2);
mean_QoS = mean(thruput_QoS+max_tpt_QoS,2);
mean_heuristic = mean(thruput_heu+max_tpt_opt,2);
mean_lin_opt = mean(thruput_lin_opt,2);
mean_slack = mean(thruput_Slack+max_tpt_slack,2);
target_rate = mean(thruput_GRS,'all');
close all;

ecdf(mean(maxrate_usrs.*rate_vals,2)); hold all;
ecdf(mean_QoS); 
%ecdf(mean_OGRS);
ecdf(mean_slack); 
%ecdf(mean_heuristic);
ecdf(mean_lin_opt);
legend('MaxRate','QoS','Slack','LinProg')
%legend('MaxRate','QoS','Deficit+mxRate','Slack','Heuristic','LinProg')

error_slack = [];
error_QoS = [];
for i = 1:n_users
    error_slack = [error_slack  sum(perform_slack(i,:)<thruput_i(i)-0.01)/n_iterations*tau_i(i)];
    error_QoS = [error_QoS   sum(perform_QoS(i,:)<thruput_i(i)-0.01)/n_iterations*tau_i(i)];
end

% ecdf(mean(maxrate_usrs.*rate_vals,2)); hold all;
% ecdf(mean_lin_opt);
% legend('MaxRate','LinProg')
toc

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
    frac_mat = reshape(x,size(ch_mat));

end

%                 while (temp_total < MAX_RBS) && (0<sum(Slack_Deficit))
%                     % Iteratively allocate RBs until no resources left or until
%                     % no service deficit
%                     temp_rate = max(rate_vals(needy_users,idx));
%                     temp_set = find(temp_rate == rate_vals(:,idx));
%                     temp_set = intersect(temp_set, needy_users);
%                     mx_id = temp_set(randi(length(temp_set)));
%     
%                     req_resources = Slack_Deficit(mx_id)/temp_rate;
%                     if req_resources > MAX_RBS
%                         M_slack(mx_id,idx) = MAX_RBS ;
%                         Slack_Deficit(mx_id) = max(Slack_Deficit(mx_id)-MAX_RBS*rate_vals(mx_id,idx),0);
%                         break;
%                     else
%                         M_slack(mx_id,idx) = req_resources;
%                         Slack_Deficit(mx_id) = 0; % No pendng Deficit 
%                         temp_total = temp_total + req_resources;
%                     end
%                 end


%rate_vals =  100*randi([10,20],n_users,n_iterations);

%----------------------------------------------------------------------------------------------
%   Recently deleted piece of code for slack algorithm
%
%             temp_rate = max(rate_vals(needy_users,idx));
%             temp_set = find(temp_rate == rate_vals(:,idx));
%             temp_set = intersect(temp_set, needy_users);
%             mx_id = temp_set(randi(length(temp_set)));
%     
%             req_resources = Slack_Deficit(mx_id)/temp_rate;
%             if req_resources < MAX_RBS
%                 % with next best user and so on
%                 M_slack(mx_id,idx) = req_resources ;
%                 Slack_Deficit(mx_id) = 0; % No pendng Deficit 
%     
%                 temp_total = req_resources;
%                 needy_users = find(Slack_Deficit);
%                 if (temp_total < MAX_RBS) && (0==sum(Slack_Deficit))
%                     idx_set = find(max_rates(idx)==rate_vals(:,idx));
%                     mx_id = idx_set(randi(length(idx_set)));
%                     M_slack(mx_id,idx) = MAX_RBS-temp_total ;
%                 else
%                     
%                     req_resources = Slack_Deficit(mx_id)/temp_rate;
%                 end
%     
%             else
%                 M_slack(mx_id,idx) = MAX_RBS ;
%                 Slack_Deficit(mx_id) = max(Slack_Deficit(mx_id)-MAX_RBS*rate_vals(mx_id,idx),0);
%             end
