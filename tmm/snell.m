function y= snell(n_1,n_2,th_1)
th_2_guess=asin(n_1*sin(th_1)/n_2);
if is_forward_angle(n_2,th_2_guess)
    y=th_2_guess;
else
    y=pi-th_2_guess;
end
end

