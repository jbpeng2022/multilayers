clc
clear
param_value = 5;

% 设置 fminbnd 的搜索范围
lower_bound = -20;
upper_bound = 10;

% 调用 fminbnd 进行优化
x_min = fminbnd(@(x) func(x, param_value), lower_bound, upper_bound)
