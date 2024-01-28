function [R,T]=coh_tmm(pol,n_list,d_list,th_0,lam_vac)
th_list=list_snell(n_list,th_0);
kz_list=2*pi*n_list.*cos(th_list)/lam_vac;
delta=kz_list.*d_list;
numlayers=length(n_list);
t_list=zeros(numlayers,numlayers);
r_list=zeros(numlayers,numlayers);
for i=1:numlayers-1
    t_list(i,i+1)=interface_t(pol, n_list(i), n_list(i+1),th_list(i), th_list(i+1));
    r_list(i,i+1)=interface_r(pol, n_list(i), n_list(i+1),th_list(i), th_list(i+1));
end
M_list=zeros(2,2,numlayers-2);
for i=2:numlayers-1
    m1=[exp(-1j*delta(i)),0;0,exp(1j*delta(i))];
    m2=[1,r_list(i,i+1);r_list(i,i+1),1];
    M_list(:,:,i-1)=1/t_list(i,i+1)*m1*m2;
end
Mtilde=[1,0;0,1];
for i=1:numlayers-2
    Mtilde=Mtilde*M_list(:,:,i);
end
m3=[1,r_list(1,2);r_list(1,2),1]/t_list(1,2);
Mtilde=m3*Mtilde;

r=Mtilde(2,1)/Mtilde(1,1);
t=1/Mtilde(1,1);

R=abs(r)^2;
T=T_from_t(pol,t,n_list(1),n_list(end),th_0,th_list(end));
end
