function R_T=tmm(thichness,nk,wavelength_all,theta)
d=thichness;%nm
R_T=zeros(length(wavelength_all),2,length(d));
d_list = [inf,d,inf];
for j=1:length(wavelength_all)
    n_list = [1,nk(j),1];
    if wavelength_all(j)/d>=0.015
        [R_T(j,1),R_T(j,2)]=coh_tmm('s',n_list, d_list, theta*pi/180, wavelength_all(j));
    else 
        c_list={'i','i','i'};
        [R_T(j,1),R_T(j,2)]=incoh_tmm('s',n_list, d_list, c_list,theta*pi/180, wavelength_all(j));
    end
end
end

