clc;
clear;
% origin file
% data_set_path='/Users/pengjiabao/Documents/multilayer/fitting/data_set/';
% material_data='Si.txt';

% the openfilters data
file_path = 'Si.txt';
fileID = fopen(file_path, 'r');
data = textscan(fileID, '%s %f %f', 'HeaderLines', 0);

openfilter_output=importdata('openfilter_output').data;

d=80000;%nm
wavelength_all=1260:1:1360;%nm
lambda_list=wavelength_all;%nm
% n_dispersion=importdata([data_set_path,material_data]).data;
% material_nk_fn = interp1(n_dispersion(:,1),n_dispersion(:,2),lambda_list,'pchip' )+1j*interp1(n_dispersion(:,1),n_dispersion(:,3),lambda_list,'pchip');

material_nk_fn = interp1(cellfun(@(x) str2num(x),data{1}),real(data{2}),lambda_list,'pchip' )+1j*interp1(cellfun(@(x) str2num(x),data{1}),abs(imag(data{2})),lambda_list,'pchip');
cellfun(@(x) str2num(x),data{1})

R_T=zeros(length(wavelength_all),2,length(d));
theta=30;
wavelength_all/d;
for i=1:length(d)
    d_list = [inf,d(i),inf];
    for j=1:length(wavelength_all)
        n_list = [1,material_nk_fn(j),1];
        if wavelength_all(j)/d>=0.015
            [R_T(j,1,i),R_T(j,2,i)]=coh_tmm('s',n_list, d_list, theta*pi/180, wavelength_all(j));
        else 
            c_list={'i','i','i'};
            [R_T(j,1,i),R_T(j,2,i)]=incoh_tmm('s',n_list, d_list, c_list,theta*pi/180, wavelength_all(j));
        end
    end
end
data_set_path='/Users/pengjiabao/Desktop/temporary data set/';
reflectance_data='Reflectance-calcs.txt';
reflectance=importdata([data_set_path,reflectance_data]).data;
transmittance_data='Transmittance-calcs.txt';
transmittance=importdata([data_set_path,transmittance_data]).data;



figure(1)
plot(wavelength_all,R_T(:,1,1),'Color','black','DisplayName','our program');
hold on 
% plot(wavelength_all,reflectance(:,2),'Color','red','LineStyle','--','DisplayName','Calculator')
% calculator

% plot(openfilter_output(:,1),openfilter_output(:,2),'Color','red','LineStyle','--','DisplayName','openfilter');


xlabel('wavelength(nm)');
ylabel('Reflection');
title(strcat('Reflection of s polarized light at',num2str(theta),'^{\circ}(black)'));
legend('show')



% figure(2)
% plot(wavelength_all,R_T(:,2,1),'Color','black','DisplayName','our program');
hold on 
% plot(wavelength_all,transmittance(:,2),'Color','red','LineStyle','--','DisplayName','Calculator')
xlabel('wavelength(nm)');
ylabel('Transmittance');
title(strcat('Transmittance of s polarized light at',num2str(theta),'^{\circ}(black)'));
legend('show')

