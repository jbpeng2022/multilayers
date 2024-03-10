function [R,T]=tmm(pol,n_list,d_list,th_0,lam_vac)
state=0;
c_list={'i'};
for i=2:length(d_list)-1
    if lam_vac/d_list(i)>=0.01
        c_list{i}='c';
    else
        c_list{i}='i';
        state=1;
    end
end
c_list{length(c_list)+1}='i';

if state
    [R,T]=incoh_tmm(pol,n_list,d_list,c_list,th_0,lam_vac);
else
    [R,T]=coh_tmm(pol,n_list,d_list,th_0,lam_vac);
end

end

