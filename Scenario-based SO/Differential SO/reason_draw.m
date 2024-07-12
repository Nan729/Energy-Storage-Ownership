clc
clear all;

beta_mat = 0.0:0.02:0.5;

load('result_diff/result_sep_disjoint')
load('result_diff/result_joint.mat')
hFig = figure;
plot(beta_mat, joint_value,'LineStyle', '-', 'LineWidth', 2, 'Color', [0.9290, 0.6940, 0.1250]	);% Plot the line
hold on; 
plot(beta_mat, value_wind,'LineStyle', '-', 'LineWidth', 2, 'Color', [0.4660, 0.6740, 0.1880]);
h_xlabel = xlabel('Differential Level \beta');
h_ylabel = ylabel('Optimization Results');
grid on;