function y=interface_r(pol,n_i,n_f,th_i,th_f)
if pol=='s'
    y=(n_i * cos(th_i) - n_f * cos(th_f)) /(n_i * cos(th_i) + n_f * cos(th_f));
elseif pol=='p'
    y=(n_f * cos(th_i) - n_i * cos(th_f)) /(n_f * cos(th_i) + n_i * cos(th_f));
end
end