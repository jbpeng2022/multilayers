function th_list=list_snell(n_list,th_0)
angles=asin(n_list(1)*sin(th_0)./n_list);
answer=is_forward_angle(n_list(1), angles(1));

if answer==1
    angles(1) = pi - angles(1);
end
answer=is_forward_angle(n_list(end), angles(end));
if answer==1
    angles(end) = pi - angles(end);
end
th_list=angles;
end