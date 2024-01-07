function [R,T]=unpolarized_RT(n_list, d_list, th_0, lam_vac)
[s_data1,s_data2] = coh_tmm('s', n_list, d_list, th_0, lam_vac);
[p_data1,p_data2] = coh_tmm('p', n_list, d_list, th_0, lam_vac);
R = (s_data1+ p_data1) / 2;
T = (s_data2 + p_data2) / 2;
% y=[R,T];
end