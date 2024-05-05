function [R,T]=tmm(pol,n_list,d_list,th_0,lam_vac)
criteria=0.01;
while 1
    state=0;
    c_list={'i'};
    for i=2:length(d_list)-1
        if lam_vac/d_list(i)>=criteria
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

    if ~isnan(R)&& ~isnan(T)
        break
    end
    criteria=criteria*2;
end

end

