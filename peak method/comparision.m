clc;
clear;
R_spectro=importdata('R_750.txt').data;
T_spectro=importdata('T_750.txt').data;



material_nk='SiO2.txt';
n_dispersion=importdata(material_nk).data;
lambda_list=R_spectro(:,1);
nk = interp1(n_dispersion(:,1),n_dispersion(:,2),lambda_list,'pchip' )+1j*interp1(n_dispersion(:,1),n_dispersion(:,3),lambda_list,'pchip');

d=fitting_fft(T_spectro,nk,lambda_list,0);