% this programm determines the thickness of films by the method shown in
% the article 《Determination of optical parameters and thickness of weakly
% absorbing thin films from reflectance and transmittance spectra》
% the air surrounding is labeled by medium 1. and the substrate is labeled
% by medium 3. 
clc;
clear;
R_spectro=importdata('R_660.txt').data;
T_spectro=importdata('T_660.txt').data;
R=R_spectro(:,2);
T=T_spectro(:,2);
lambdalist=R_spectro(:,1);
theta=0;
n_1=1;
n_3=1;
%% first stage
% step1 get (R_,T+,lambda) and (R+,T_,lambda) from exp
[~,locs_R_min] = findpeaks(-R);
[~,locs_R_max] = findpeaks(R);
[~,locs_T_min] = findpeaks(-T);
[~,locs_T_max] = findpeaks(T);
R_T_1=[R(locs_R_max),T(locs_T_min),lambdalist(locs_R_max)];%(R_,T+,lambda)
R_T_2=[R(locs_R_min),T(locs_T_max),lambdalist(locs_R_min)];

% first calculation before lteration
eta_1=n_1*cos(theta);
eta_3=n_3*cos(theta);
eta_2_imag=0;
A_posi=(1+sqrt(R_T_1(:,1)))./sqrt(R_T_1(:,2));
B_posi=(1-sqrt(R_T_1(:,1)))./sqrt(R_T_1(:,2));
eta_2_real=sqrt(eta_1*eta_3*(eta_3+eta_1*A_posi.^2)./(eta_1+eta_3*B_posi.^2));
h_before=-100000;

while 1
% calculate delta_2
    eta_2=eta_2_real-1j*eta_2_imag;
    delta_2=angle((eta_2-eta_3)./(eta_2+eta_3));
% step 4 calculate h
    pairs=0;
    h_lk=0;
    for l=1:length(delta_2)-1
        for k=l+1:length(delta_2)
            p=k-l;
            lambda_l=R_T_1(l,3);
            lambda_k=R_T_1(k,3);
            eta_2_real_k=eta_2_real(k);
            eta_2_real_l=eta_2_real(l);
            h_lk=h_lk+(2*p*pi+delta_2(l)-delta_2(k))/4/pi*lambda_l*lambda_k/(eta_2_real_k*lambda_l-eta_2_real_l*lambda_k);
            pairs=pairs+1;
        end
    end
    h=abs(h_lk/pairs)
% step 5 calculate alpha
    C=(1-R_T_1(:,1))./R_T_1(:,2);
    alpha1=2*abs(eta_2).^2*eta_3.*C+sqrt((2*abs(eta_2).^2*eta_3.*C).^2+eta_2_real.^2.*((abs(eta_2).^2+eta_3^2).^2-(2*eta_2_real*eta_3).^2));
    alpha2=eta_2_real.*((eta_2_real+eta_3).^2+eta_2_imag.^2);
    alpha=log(alpha1./alpha2);
    k_0=2*pi./R_T_1(:,3);
    eta_2_imag=alpha/2./k_0/h;
% R_b,R,T
    eta_2=eta_2_real-1j*eta_2_imag;
    r1=abs((eta_1-eta_2)./(eta_1+eta_2));
    t1=abs(2*eta_1./(eta_1+eta_2));
    r2=abs((eta_2-eta_3)./(eta_2+eta_3));
    t2=abs(2*eta_2./(eta_2+eta_3));
    r3=abs((eta_3-eta_1)/(eta_3+eta_1));
    t3=abs(2*eta_3/(eta_1+eta_3));

    R_b=abs((r2+r1.*exp(-2j*k_0*h.*eta_2))./(1+r2.*r1.*exp(-2j*k_0*h.*eta_2))).^2;
    R_T_1(:,1)=R_T_1(:,1)-r3^2*R_T_1(:,2).^2./(1-r3^2.*R_b);
    R_T_1(:,2)=eta_3*(1-r3^2*R_b)/eta_1/t3^2.*R_T_1(:,2);
    % R=R-r3^2*exp(-2*k_0*h_3*eta_3_imag)*T^2/(1-r3^2*exp(-2*k_0*h_3*eta_3_imag)*R_b);
    % T=eta_3*(1-r3^2*exp(-2*k_0*h_3*eta_3_imag)*R_b)/eta_1/t3^2/exp(-k_0*h_3*eta_3_imag)*T;    
% eta_2_real
    eta_2_real=eta_3*(1-abs(r1).^2).*abs(t2).^2.*exp(-alpha)./R_T_1(:,2)./(1+r1.*r2.*exp(-alpha)).^2;
    if abs(h_before-h)<=5             %nm
        break
    end
    h_before=h;
end
eta_2=eta_2_real-1j*eta_2_imag;
%% second stage
delta_lambda=1;
con_all=abs([eta_2_real;0]-[0;eta_2_real]);
con_all(1)=con_all(2);
con_all(end)=con_all(end-1);
n_dispersion=[];
for i=1:length(R(locs_R_max))
   eta_2_1=eta_2(i);
   eta_2_real_0=real(eta_2_1);
   n_dispersion=[n_dispersion;R_T_1(i,3),real(eta_2_1)/cos(theta),imag(eta_2_1)/cos(theta)];
   if  i==1 %first peaks left direction
       con=con_all(1);
       lambda=R_T_1(i,3)-delta_lambda;
       eta_2_1_imag=imag(eta_2_1);
       eta_2_1_real=real(eta_2_1);
       C=(1-R_spectro(R_spectro(:,1)==lambda,2))/T_spectro(T_spectro(:,1)==lambda,2);
       alpha1=2*abs(eta_2_1).^2*eta_3.*C+sqrt((2*abs(eta_2_1).^2*eta_3.*C).^2+eta_2_1_real.^2.*((abs(eta_2_1).^2+eta_3^2).^2-(2*eta_2_1_real*eta_3).^2));
       alpha2=eta_2_1_real*((eta_2_1_real+eta_3)^2+eta_2_1_imag^2);
       alpha=log(alpha1/alpha2);
       while 1
            R_mea=R_spectro(R_spectro(:,1)==lambda,2);
            for k=1:2
                r1=abs((eta_1-eta_2_1)/(eta_1+eta_2_1));
                r2=abs((eta_2_1-eta_3)/(eta_2_1+eta_3));
                delta=angle((eta_2_1-eta_3)/(eta_2_1+eta_3));
                k_0=2*pi/lambda;
                fun=@(x)((r1^2+r2^2*exp(-2*alpha)+2*r1*r2*exp(-alpha)*cos(2*k_0*h*x+delta))/(1+r1^2*r2^2*exp(-2*alpha)+2*r1*r2*exp(-alpha)*cos(2*k_0*h*x+delta))-R_mea)^2;
                x0 =eta_2_real_0;
                x= fminsearchbnd(fun,x0,x0-con,x0+con);
                eta_2_1_real=x;
                eta_2_1=eta_2_1_real-1j*eta_2_1_imag;
                C=(1-R_spectro(R_spectro(:,1)==lambda,2))/T_spectro(T_spectro(:,1)==lambda,2);
                alpha1=2*abs(eta_2_1).^2*eta_3.*C+sqrt((2*abs(eta_2_1).^2*eta_3.*C).^2+eta_2_1_real.^2.*((abs(eta_2_1).^2+eta_3^2).^2-(2*eta_2_1_real*eta_3).^2));
                alpha2=eta_2_1_real.*((eta_2_1_real+eta_3).^2+eta_2_1_imag.^2);
                alpha=log(alpha1/alpha2);
                eta_2_1_imag=alpha/2/k_0/h;
                eta_2_1=eta_2_1_real-1j*eta_2_1_imag;
                % eta_2_real_0=real(eta_2_1);
            end
            n_dispersion=[n_dispersion;lambda,eta_2_1_real/cos(theta),eta_2_1_imag/cos(theta)];
            if lambda==R_spectro(1,1)
                break
            end
            lambda=lambda-delta_lambda;
       end
       %right direction
       con=con_all(i+1);
       lambda=R_T_1(i,3)+delta_lambda;
       eta_2_1_imag=imag(eta_2_1);
       eta_2_1_real=real(eta_2_1);
       C=(1-R_spectro(R_spectro(:,1)==lambda,2))/T_spectro(T_spectro(:,1)==lambda,2);
       alpha1=2*abs(eta_2_1).^2*eta_3.*C+sqrt((2*abs(eta_2_1).^2*eta_3.*C).^2+eta_2_1_real.^2.*((abs(eta_2_1).^2+eta_3^2).^2-(2*eta_2_1_real*eta_3).^2));
       alpha2=eta_2_1_real.*((eta_2_1_real+eta_3).^2+eta_2_1_imag.^2);
       alpha=log(alpha1/alpha2);
       while 1 
            R_mea=R_spectro(R_spectro(:,1)==lambda,2);
            for k=1:2
                r1=abs((eta_1-eta_2_1)/(eta_1+eta_2_1));
                r2=abs((eta_2_1-eta_3)/(eta_2_1+eta_3));
                delta=angle((eta_2_1-eta_3)/(eta_2_1+eta_3));
                k_0=2*pi/lambda;
                fun=@(x)((r1^2+r2^2*exp(-2*alpha)+2*r1*r2*exp(-alpha)*cos(2*k_0*h*x+delta))/(1+r1^2*r2^2*exp(-2*alpha)+2*r1*r2*exp(-alpha)*cos(2*k_0*h*x+delta))-R_mea)^2;
                x0 =eta_2_real_0;
                x= fminsearchbnd(fun,x0,x0-con,x0+con);
                eta_2_1_real=x;
                eta_2_1=eta_2_1_real-1j*eta_2_1_imag;
                C=(1-R_spectro(R_spectro(:,1)==lambda,2))/T_spectro(T_spectro(:,1)==lambda,2);
                alpha1=2*abs(eta_2_1).^2*eta_3.*C+sqrt((2*abs(eta_2_1).^2*eta_3.*C).^2+eta_2_1_real.^2.*((abs(eta_2_1).^2+eta_3^2).^2-(2*eta_2_1_real*eta_3).^2));
                alpha2=eta_2_1_real.*((eta_2_1_real+eta_3).^2+eta_2_1_imag.^2);
                alpha=log(alpha1/alpha2);
                eta_2_1_imag=alpha/2/k_0/h;
                eta_2_1=eta_2_1_real-1j*eta_2_1_imag;
                % eta_2_real_0=real(eta_2_1);
            end
            n_dispersion=[n_dispersion;lambda,real(eta_2_1)/cos(theta),imag(eta_2_1)/cos(theta)];
            if lambda==R_T_1(2,3)
                break
            end
            lambda=lambda+delta_lambda;
       end
   else%right direction
       con=con_all(i+1);
       lambda=R_T_1(i,3)+delta_lambda;
       eta_2_1_imag=imag(eta_2_1);
       eta_2_1_real=real(eta_2_1);
       C=(1-R_spectro(R_spectro(:,1)==lambda,2))/T_spectro(T_spectro(:,1)==lambda,2);
       alpha1=2*abs(eta_2_1).^2*eta_3.*C+sqrt((2*abs(eta_2_1).^2*eta_3.*C).^2+eta_2_1_real.^2.*((abs(eta_2_1).^2+eta_3^2).^2-(2*eta_2_1_real*eta_3).^2));
       alpha2=eta_2_1_real.*((eta_2_1_real+eta_3).^2+eta_2_1_imag.^2);
       alpha=log(alpha1/alpha2);
       while 1 
            R_mea=R_spectro(R_spectro(:,1)==lambda,2);
            for k=1:3
                r1=abs((eta_1-eta_2_1)/(eta_1+eta_2_1));
                r2=abs((eta_2_1-eta_3)/(eta_2_1+eta_3));
                delta=angle((eta_2_1-eta_3)/(eta_2_1+eta_3));
                k_0=2*pi/lambda;
                fun=@(x)((r1^2+r2^2*exp(-2*alpha)+2*r1*r2*exp(-alpha)*cos(2*k_0*h*x+delta))/(1+r1^2*r2^2*exp(-2*alpha)+2*r1*r2*exp(-alpha)*cos(2*k_0*h*x+delta))-R_mea)^2;
                x0=eta_2_real_0;
                x= fminsearchbnd(fun,x0,x0-con,x0+con);
                eta_2_1_real=x;
                eta_2_1=eta_2_1_real-1j*eta_2_1_imag;
                C=(1-R_spectro(R_spectro(:,1)==lambda,2))/T_spectro(T_spectro(:,1)==lambda,2);
                alpha1=2*abs(eta_2_1).^2*eta_3.*C+sqrt((2*abs(eta_2_1).^2*eta_3.*C).^2+eta_2_1_real.^2.*((abs(eta_2_1).^2+eta_3^2).^2-(2*eta_2_1_real*eta_3).^2));
                alpha2=eta_2_1_real.*((eta_2_1_real+eta_3).^2+eta_2_1_imag.^2);
                alpha=log(alpha1/alpha2);
                eta_2_1_imag=alpha/2/k_0/h;
                eta_2_1=eta_2_1_real-1j*eta_2_1_imag;
            end
            n_dispersion=[n_dispersion;lambda,real(eta_2_1)/cos(theta),imag(eta_2_1)/cos(theta)];
            
            if i ~=length(R_T_1(:,1))
                if lambda==R_T_1(i+1,1)
                    break
                end
            end

            if lambda==R_spectro(end,1)
                break
            end
            lambda=lambda+delta_lambda;
       end
   end
end
%% plot the results
n_dispersion=sortrows(n_dispersion,1);
% p= polyfit(n_dispersion(:,1),n_dispersion(:,2),6);
% n_dispersion_poly= polyval(p,n_dispersion(:,1));
f = fit(n_dispersion(:,1),n_dispersion(:,2),'a+b/x^2+c/x^4');


% scatter(n_dispersion(:,1),n_dispersion(:,2))
% plot(n_dispersion(:,1),n_dispersion(:,2))
% hold on 
% plot(n_dispersion(:,1),n_dispersion_poly,'LineStyle','--','LineWidth',3)
plot( f,n_dispersion(:,1) , n_dispersion(:,2)) 
hold on 
% plot(n_dispersion(:,1),n_dispersion(:,3),'LineWidth',3)
hold on 

material=importdata('SiO2.txt').data;
plot(material(:,1),material(:,2),'LineWidth',2)






