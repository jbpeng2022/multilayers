function merit= spectro_cal(eta_2_real,R_mea,k_0,h,eta_1,eta_3,eta_2_imag_0,C)
eta_2=eta_2_real-1j*eta_2_imag_0;
r1=abs((eta_1-eta_2)/(eta_1+eta_2));
r2=abs((eta_2-eta_3)/(eta_2+eta_3));
delta2=angle((eta_2-eta_3)/(eta_2+eta_3));
delta1=angle((eta_1-eta_2)/(eta_1+eta_2));
% delta1=pi;

alpha1=2*abs(eta_2).^2*eta_3.*C+sqrt((2*abs(eta_2).^2*eta_3.*C).^2+eta_2_real.^2.*((abs(eta_2).^2+eta_3^2).^2-(2*eta_2_real*eta_3).^2));
alpha2=eta_2_real.*((eta_2_real+eta_3).^2+eta_2_imag_0.^2);
alpha=log(alpha1/alpha2);

R_cal=(r1^2+r2^2*exp(-2*alpha)+2*r1*r2*exp(-alpha)*cos(2*k_0*h*eta_2_real+delta2-delta1))/(1+r1^2*r2^2*exp(-2*alpha)+2*r1*r2*exp(-alpha)*cos(2*k_0*h*eta_2_real+delta2+delta1));
merit=(R_cal-R_mea)^2;
end
