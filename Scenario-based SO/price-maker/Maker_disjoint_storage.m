function [outputArg1,outputArg2, x_D] = Maker_disjoint_storage(w_R)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
DA_path = '../../inter_data/DA_price_Aug_HB_SOUTH.mat';
w_path = '../../inter_data/wind_fore_summer.mat';
save_path = './result/disjoint_SO_maker.mat';
tempname = '../inter/wind_data.mat';

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

% use external training set
external_RT_name = '../inter/RT_gap_train.mat';
p_R = (1+load(external_RT_name).data) .* p;
sample_num = size(p_R);
S = sample_num(2);

load(tempname, 'z');

slope = 0.8/1000;

% storage model
cvx_begin quiet
    variable x_D(T);
    variable sta_D(T);
    variable x_R(T,S);
    variable dis_R(T,S);
    variable sta_R(T,S);
    expression x_obj;

    x_obj = p' * (x_D);

    for s = 1:S
        x_obj = x_obj + 1/S* ( ( p_R(:, s) - slope * w_R(:, s)  )' * x_R(:, s) - slope * x_R(:, s)' * x_R(:, s)  );
    end
    
    % Objective function
    maximize(x_obj)

    % Constraints
    subject to
        sta_D == add_mat * sta_D - (1/eta_d)*x_D + s_init
        sta_D <= s_ub
        sta_D >= s_lb
        x_D >= -ramp
        x_D <= ramp

        for s = 1:S
            x_D+ x_R(:, s) == dis_R(:, s)
            sta_R(:, s) == add_mat * sta_R(:, s) - (1/eta_d)*dis_R(:, s) + s_init
            sta_R(:, s) <= s_ub
            sta_R(:, s) >= s_lb
            dis_R(:, s) >= -ramp
            dis_R(:, s) <= ramp
        end
cvx_end

% Record the optimal value
optimal_value_storage = cvx_optval;


outputArg1 = optimal_value_storage;
outputArg2 = x_R;
end