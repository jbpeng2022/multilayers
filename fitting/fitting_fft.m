function [d_fft,d_interval]=fitting_fft(R_v,nk,wavelength_all,theta)
%% fft
% N=128;% the number of sampling
v_min=1/R_v(length(R_v),1);%nm
v_max=1/R_v(1,1);
% R_t=fftshift(fft(R_v(:,2).*hamming(length(R_v(:,2)),'periodic'),N));
% R_t=fftshift(fft(R_v(:,2),N));
% plot(-N/2:N/2-1,abs(R_t).^2);

% N=101; %this should pay more attention! how to choose the suitable points is important
N=length(R_v(:,2));
% N=length(R_v(:,2));
R_t=fftshift(fft(R_v(:,2).*hamming(length(R_v(:,2)),'periodic'),N));
%% search for maximum magnitude in the power spectrum
index=-N/2:N/2-1;
[pks,locs] = findpeaks(abs(R_t).^2);
sortedR_t = sort(pks, 'descend');
m=index(locs(pks==sortedR_t(2)));
m=sum(abs(m))/2;
%% determine the thickness
d_fft=m/2/(real(nk(1))*v_max-real(nk(length(nk)))*v_min);
%% applicable analysis
d_min=1/2/(real(nk(1))*v_max-real(nk(length(nk)))*v_min);
d_max=length(R_t)/2*d_min;
if d_fft>d_max || d_fft<d_min
    disp('the prediction may be inaccurate')
end
%% grid search 
M=32;
d_thickness_prediction=linspace(d_fft-3*d_min,d_fft+3*d_min,M);
d_interval=6*d_min/(M-1);
R_T_all=cell(M,1);
for i=1:M
   d=d_thickness_prediction(i);
   R_T=tmm(d,nk,wavelength_all,theta);
   R_T_all{i}=R_T;
end
% calculate all mse then find the minia
chi_error=cellfun(@(x) sum((x(:,1)-R_v(:,2)).^2),R_T_all);
[~,I]=min(chi_error);
d_fft=d_thickness_prediction(I);
end