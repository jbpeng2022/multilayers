% the program is going to take a experiment to verify the applicable of
% ffting program combined with fft and levenberg marqudat algorithm.
clc;
clear;
%% import data from txt
data_set_path='/Users/pengjiabao/Documents/multilayer/fitting/data_set/';
target_data='Reflectance_1_4_1.txt';
R_v=importdata([data_set_path,target_data]).data;
material_data='Si.txt';
lambda_list=R_v(:,1);%nm
n_dispersion=importdata([data_set_path,material_data]).data;
material_nk_fn = interp1(n_dispersion(:,1),n_dispersion(:,2),lambda_list,'pchip' )+1j*interp1(n_dispersion(:,1),n_dispersion(:,3),lambda_list,'pchip');
%% apply fft to determine the thickness
theta=0;
[d_fft,d_interval]=fitting_fft(R_v,material_nk_fn,lambda_list,theta);
optimized_thickness=fitting_regresion_analysis_lm(d_fft,d_interval,material_nk_fn,lambda_list,theta,R_v(:,2))