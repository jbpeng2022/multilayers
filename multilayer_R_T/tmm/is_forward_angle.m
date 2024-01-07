function answer=is_forward_angle(n,theta)
ncostheta=n*cos(theta);
epsilon=2.220446049250313e-16;
if abs(ncostheta)>100*epsilon
    answer=imag(ncostheta)>0;
else
    answer=real(ncostheta)>0;
end

end