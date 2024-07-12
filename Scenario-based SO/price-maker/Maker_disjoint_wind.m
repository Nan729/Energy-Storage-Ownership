function [outputArg1,outputArg2, w_D] = Maker_disjoint_wind(x_R)
%UNTITLED2 此处显示有关此函数的摘要
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

% disjoint wind original
cvx_begin quiet
    variable w_D(T);
    variable w_R(T,S);
    
    expression w_obj
    w_obj = p' * (w_D);

    for s = 1:S
        w_obj = w_obj + 1/S* ( ( p_R(:, s) - slope * x_R(:, s) )' * w_R(:, s) - slope * w_R(:, s)' * w_R(:, s)  );
    end
    
    % Objective function
    maximize(w_obj)

    % Constraints
    subject to
        w_D >= 0
        w_D <= w_cap

        for s = 1:S
            w_D+ w_R(:, s) == z(:, s)
        end
cvx_end

optimal_value_wind = cvx_optval;




outputArg1 = optimal_value_wind;
outputArg2 = w_R;
end