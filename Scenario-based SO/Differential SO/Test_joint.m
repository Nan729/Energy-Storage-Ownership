function [outputArg1,outputArg2] = Test_joint(beta, w_D, x_D)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
DA_path = '../../inter_data/DA_price_Aug_HB_SOUTH.mat';
tempname = '../inter/wind_data_test.mat';
w_path = '../../inter_data/wind_fore_summer.mat';

data = load(DA_path);
p = data.data';
T = length(p);

% wind output
w_cap = 1000;
wind_min = load(w_path);

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

cvx_quiet(true);

value = zeros(S, 1);
for s = 1:S
    cvx_begin
        variable o_R(T);
        variable dis_R(T);
        variable sta_R(T);
        variable profit(1);

        % Objective function
        maximize(profit)

        subject to
            w_D+ x_D+ o_R== dis_R + z(:, s)
            sta_R == add_mat * sta_R - (1/eta_d)*dis_R + s_init
            sta_R <= s_ub
            sta_R >= s_lb
            dis_R >= -ramp
            dis_R <= ramp

            profit <= p_R_sell(:,s)'*o_R
            profit <= p_R_buy(:,s)'*o_R

    cvx_end

    value(s) = cvx_optval;
end

value_new = p' * (x_D+w_D) + value;

avg = sum(value_new)/S;
std_j = std(value_new);

outputArg1 = avg/1000;
outputArg2 = std_j/1000;
end