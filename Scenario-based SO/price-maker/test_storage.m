function [outputArg1,outputArg2] = test_storage(s, w_R)
%UNTITLED7 此处显示有关此函数的摘要
%   此处显示详细说明
DA_path = '../../inter_data/DA_price_Aug_HB_SOUTH.mat';
save_path = './result/disjoint_SO_maker.mat';
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

load(tempname, 'z');

load(save_path, "w_D", "x_D");

slope = 0.8/1000;

% Suppress intermediate output
cvx_quiet(true);

cvx_begin
    variable x_R(T);
    variable dis_R(T);
    variable sta_R(T);

    % Objective function
    maximize(p' * (x_D) + (p_R(:,s) - slope * w_R)' * x_R - slope * x_R' * x_R)

    % Constraints
    subject to
        x_D+ x_R == dis_R
        sta_R == add_mat * sta_R - (1/eta_d)*dis_R + s_init
        sta_R <= s_ub
        sta_R >= s_lb
        dis_R >= -ramp
        dis_R <= ramp
cvx_end

value_storage = cvx_optval;


outputArg1 = value_storage;
outputArg2 = x_R;
end