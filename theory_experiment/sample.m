clc;
clear;
%% low doping concentration
% where N is the doping concentration
% d is the thickness of single layer / should be normalized to unit with nm
% theta is the incident angle. 
% the unit of wavelength should be m and nm in doped_silicon function and
% coh_tmm function, respectively.
data_set_path='/Users/pengjiabao/Documents/multilayer/fitting/data_set/';
material_data='Si.txt';

N=1e19;%cm-3
% d=linspace(0.1,10,5)*10^3;%nm
d=1000;%nm
wavelength_all=200:1:1000;%nm
lambda_list=wavelength_all;%nm
n_dispersion=importdata([data_set_path,material_data]).data;
material_nk_fn = interp1(n_dispersion(:,1),n_dispersion(:,2),lambda_list,'pchip' )+1j*interp1(n_dispersion(:,1),n_dispersion(:,3),lambda_list,'pchip');

R_T=zeros(length(wavelength_all),2,length(d));
theta=0;
for i=1:length(d)
    d_list = [inf,d(i),inf];
    for j=1:length(wavelength_all)
        % [n,k]=doped_silicon(N,'n',wavelength_all(j)*10^-9);
        % n_list = [1,n+1i*k,1];
        n_list = [1,material_nk_fn(j),1];
        [R_T(j,1,i),R_T(j,2,i)]=coh_tmm('s',n_list, d_list, theta*pi/180, wavelength_all(j));
    end
end
for j=1:length(d)
    figure(1)
    plot(wavelength_all,R_T(:,1,j),'Color','black');
    hold on 
    xlabel('wavelength(nm)');
    ylabel('Fraction reflected');
    title(strcat('Reflection of s polarized light at',num2str(theta),'^{\circ}(black)'));
    figure(2)
    plot(wavelength_all,R_T(:,2,j),'Color','black');
    hold on 
    xlabel('wavelength(nm)');
    ylabel('Fraction of power transmitted');
    title(strcat('Rafraction of s polarized light at',num2str(theta),'^{\circ}(black)'));
end
%% high doping concentration
% N=1*10^21;
% d=400*10^3;%nm
% wavelength_all=linspace(200,2000,10000);%nm
% R_T=zeros(length(wavelength_all),2);
% theta=45;
% d_list = [inf,d,inf];
% for j=1:length(wavelength_all)
%     [n,k]=doped_silicon(N,'p',wavelength_all(j)*10^-9);
%     n_list=[1,n+1i*k,1];
%     [R_T(j,1),R_T(j,2)]=coh_tmm('s',n_list, d_list, theta*pi/180, wavelength_all(j));
% end
% figure(3)
% plot(wavelength_all,R_T(:,1),'Color','black');
% hold on 
% xlabel('wavelength(nm)');
% ylabel('Fraction reflected');
% title(strcat('Reflection of s polarized light at',num2str(theta),'^{\circ}(black)'));
% figure(4)
% plot(wavelength_all,R_T(:,2),'Color','black');
% hold on 
% xlabel('wavelength(nm)');
% ylabel('Fraction of power transmitted');
% title(strcat('Rafraction of s polarized light at',num2str(theta),'^{\circ}(black)'));





