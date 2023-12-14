function y=interface_t(pol, n_i, n_f, th_i, th_f)
if pol == 's'
    y=2 * n_i * cos(th_i) / (n_i * cos(th_i) + n_f * cos(th_f));
elseif pol == 'p'
    y=2 * n_i * cos(th_i) / (n_f * cos(th_i) + n_i * cos(th_f));
end