% where N is the doped_concentration, and type is the dope_type 
% and T is the temperture in kelvin, omega is radius frequency of light
function [n,k]=doped_silicon(N,type,wavelength)
%% calculate the carrier concentration
T=293;%in kelvin
v=2.99792458*10^8;%light speed vacuum
theta=T/300;
omega=2*pi*v/wavelength;
if type=='n'
    A=0.0824*theta^-1.622;
    N_0=1.6*10^18*theta^0.7267;
    if N<N_0
        B=0.4722*theta^0.0652;
    elseif N>=N_0
        B=1.23-0.3162*theta;
    end
elseif type=='p'
    A=0.2364*theta^-1.474;
    N_0=1.577*10^18*theta^0.46;
    if N<N_0
        B=0.433*theta^0.2213;
    elseif N>=N_0
        B=1.268-0.338*theta;
    end
else 
    print('please input type with n or p')
    quit
end
degree_of_ionization=1-A*exp(-(B*log(N/N_0))^2);
carrier_concentration=N*degree_of_ionization;
%% calculate the mobility
if type=='n'
    mu_1=68.5;
    mu_max=1414;
    mu_2=56.1;
    C_r=9.2*10^16;
    C_s=3.41*10^20;
    alpha=0.711;
    beta=198;
    N_e=carrier_concentration;
    mu=mu_1+(mu_max-mu_1)/(1+(N_e/C_r)^alpha)-mu_2/(1+(C_s/N_e)^beta);
elseif type=='p'
    mu_1=44.9;
    mu_max=470.5;
    mu_2=29;
    C_r=2.23*10^17;
    C_s=6.10*10^20;
    alpha=0.719;
    beta=2;
    p_c=9.23*10^16;
    N_h=carrier_concentration;
    mu=mu_1*exp(-p_c/N_h)+mu_max/(1+(N_h/C_r)^alpha)-mu_2/(1+(C_s/N_h)^beta);
else 
    print('please input type with n or p')
    quit
end
%% by using the values for carrier concentration and mobility, calculate the optical constant
e=1.602*10^-19;%-10esu
m_0=9.109*10^-31;%KG
epsilion0=8.854*10^-12;%ui C^2/N M^2
dielectric_limitation=11.78;
if type=='n'
    relative_mass=0.27*m_0;
elseif type=='p'
    relative_mass=0.37*m_0;
else
    print('please input type with n or p')
    quit
end
omega_p=sqrt(carrier_concentration*e^2/relative_mass/epsilion0)*10^3;
gamma=e/relative_mass/mu*10^4;
relative_permittivity=dielectric_limitation-omega_p^2/omega/(omega+1i*gamma);
syms u v
eqns = [u^2-v^2 == real(relative_permittivity), 2*u*v== imag(relative_permittivity)];
S = solve(eqns,[u v]);
U=double(S.u);
V=double(S.v);
n=U(U>0);
k=V(V>0);
end








