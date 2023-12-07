function T=T_from_t(pol,t,n_i,n_f,th_i,th_f)
if pol=='s'
    T=abs(t^2)*(real(n_f*cos(th_f))/real(n_i*cos(th_i)));
elseif pol=='p'
    T=abs(t^2)*(real(n_f*conj(cos(th_f)))/real(n_i*conj(cos(th_i))));
end