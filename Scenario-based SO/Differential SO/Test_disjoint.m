function [outputArg1,outputArg2] = Test_disjoint(beta, w_D, x_D)
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明
DA_path = '../../inter_data/DA_price_Aug_HB_SOUTH.mat';
tempname = '../inter/wind_data_test.mat';


data = load(DA_path);
p = data.data';
T = length(p);
% wind output
w_cap = 1000;

s_cap = 1000;
ramp_rate = 0.2;
ramp = ramp_rate*s_cap;
s_ub = ones(T,1)*s_cap;
s_ub(T) = s_ub(T) - 0.5*s_cap;
s_lb = zeros(T,1);
s_lb(T) = s_lb(T) + 0.5*s_cap;
s_init = zeros(T,1);
s_init(1) = 0.5*s_cap;
eta_d = 1;
add_mat = tril(ones(T), -1)- tril(ones(T), -2);

% use external testing set
external_RT_name = '../inter/RT_gap_test.mat';
p_R = (1+load(external_RT_name).data) .* p;
sample_num = size(p_R);
S = sample_num(2);

p_R_sell = p_R * (1-beta);
p_R_buy = p_R * (1+beta);

load(tempname, 'z');

% Suppress intermediate output
cvx_quiet(true);

% disjoint wind
value_wind = zeros(S,1);
for s = 1:S
    cvx_begin
        variable w_R(T);
        variable profit(1);
        % Objective function
        maximize(p' * w_D + profit)

        % Constraints
        subject to
            w_D+ w_R == z(:, s)
            
            profit <= (p_R_sell(:,s))' * w_R
            profit <= (p_R_buy(:,s))' * w_R
    cvx_end
    value_wind(s) = cvx_optval;
end

avg_wind = mean(value_wind)/1000;
std_wind = std(value_wind)/1000;

% storage model
value_storage = zeros(S,1);
for s = 1:S
    cvx_begin
        variable x_R(T);
        variable dis_R(T);
        variable sta_R(T);
        variable profit_2(1);

        % Objective function
        maximize(p' * (x_D) + profit_2)

        % Constraints
        subject to
            x_D+ x_R == dis_R
            sta_R == add_mat * sta_R - (1/eta_d)*dis_R + s_init
            sta_R <= s_ub
            sta_R >= s_lb
            dis_R >= -ramp
            dis_R <= ramp

            profit_2 <= (p_R_sell(:,s))' * x_R
            profit_2 <= (p_R_buy(:,s))' * x_R
    cvx_end

    value_storage(s) = cvx_optval;
end

avg_storage = mean(value_storage)/1000;
std_storage = std(value_storage)/1000;

avg_sum = mean(value_storage+value_wind)/1000;
std_sum = std(value_storage+value_wind)/1000;


outputArg1 = avg_sum;
outputArg2 = std_sum;

end