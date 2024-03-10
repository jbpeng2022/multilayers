% the function is going to calculate the case of incoherent which is a more
% common case.
function [R,T]= incoh_tmm(pol,n_list,d_list,c_list,th_0,lam_vac)
%% stacks division 
inc_from_all=zeros(length(n_list),1); 
stack_from_all=cell(length(n_list),1);
inc_index=1;
stack_index=1;
all_from_stack={};
all_from_inc=[];
stack_d_list={};
stack_n_list={};
stack_from_inc=[];

stack_in_progress=false;
for alllayer_index=1:length(n_list)
    if c_list{alllayer_index}=='c' %this layer is coherent.
        inc_from_all(alllayer_index)=nan;
        if ~stack_in_progress %this layer is starting new stack
            stack_in_progress=true;
            ongoing_stack_d_list=[inf,d_list(alllayer_index)];
            ongoing_stack_n_list=[n_list(alllayer_index-1),n_list(alllayer_index)];
            stack_from_all{alllayer_index}=[stack_index,1];
            all_from_stack{stack_index}=alllayer_index;%%
            within_stack_index=1;
        else
            ongoing_stack_d_list=[ongoing_stack_d_list,d_list(alllayer_index)];
            ongoing_stack_n_list=[ongoing_stack_n_list,n_list(alllayer_index)];
            within_stack_index=within_stack_index+1;
            stack_from_all{alllayer_index}=[stack_index,within_stack_index];
            all_from_stack{stack_index}=[all_from_stack{stack_index},alllayer_index];%%
        end
    elseif c_list{alllayer_index}=='i'
        stack_from_all{alllayer_index}=nan;
        inc_from_all(alllayer_index)=inc_index;
        all_from_inc=[all_from_inc,alllayer_index];
        if ~stack_in_progress %previous layer is also incoherent ！！
            stack_from_inc=[stack_from_inc,nan];
        else %previous layer is coherent
            stack_in_progress=false;
            stack_from_inc=[stack_from_inc,stack_index];
            ongoing_stack_d_list=[ongoing_stack_d_list,inf];
            stack_d_list{stack_index}=ongoing_stack_d_list;
            ongoing_stack_n_list=[ongoing_stack_n_list,n_list(alllayer_index)];
            stack_n_list{stack_index}=ongoing_stack_n_list;
            all_from_stack{stack_index}=[all_from_stack{stack_index},alllayer_index];%%
            stack_index=stack_index+1;
        end
        inc_index=inc_index+1;
    end
end
num_inc_layers = length(all_from_inc);
num_stacks = length(all_from_stack);

%% caculate the R and T by using the transfer matrix
th_list = list_snell(n_list, th_0);
coh_tmm_data_list =zeros(num_stacks,2);
coh_tmm_bdata_list =zeros(num_stacks,2);
for i=1:num_stacks
    coh_tmm_data_list(i,:)=coh_tmm(pol,stack_n_list{i},stack_d_list{i},th_list(all_from_stack{i}(1)),lam_vac);%the format of all_from_stack should be checked again
    th_f = snell(n_list(1), n_list(end), th_0);
    coh_tmm_bdata_list(i,:)=coh_tmm(pol,flip(stack_n_list{i}),flip(stack_d_list{i}),th_f,lam_vac);
end
P_list=zeros(num_inc_layers);
for inc_index=2:num_inc_layers-1 %skip the fist and the last entries of inf
    i=all_from_inc(inc_index);
    P_list(inc_index)=exp(-4*pi*d_list(i)*imag((n_list(i)* cos(th_list(i)))/lam_vac));
    if P_list(inc_index)<1e-30
        P_list(inc_index)=1e-30;
    end
end
T_list=zeros(num_inc_layers,num_inc_layers);
R_list=zeros(num_inc_layers,num_inc_layers);
for inc_index=1:num_inc_layers-1
    alllayer_index=all_from_inc(inc_index);
    nextstack_index=stack_from_inc(inc_index+1);
    if isnan(nextstack_index) %next layer is incoherent
        R_list(inc_index,inc_index+1)=abs(interface_r(pol,n_list(alllayer_index),n_list(alllayer_index+1),th_list(alllayer_index),th_list(alllayer_index+1)))^2;
        T_list(inc_index, inc_index+1)=abs(interface_t(pol,n_list(alllayer_index),n_list(alllayer_index+1),th_list(alllayer_index),th_list(alllayer_index+1)))^2;
        R_list(inc_index+1, inc_index)= abs(interface_r(pol, n_list(alllayer_index+1),n_list(alllayer_index),th_list(alllayer_index+1),th_list(alllayer_index)))^2;
        T_list(inc_index+1, inc_index)=abs(interface_t(pol, n_list(alllayer_index+1),n_list(alllayer_index),th_list(alllayer_index+1),th_list(alllayer_index)))^2;
    else %next layer is coherent
        R_list(inc_index,inc_index+1)=coh_tmm_data_list(nextstack_index,1);
        T_list(inc_index,inc_index+1)=coh_tmm_data_list(nextstack_index,2);
        R_list(inc_index+1,inc_index)=coh_tmm_bdata_list(nextstack_index,1);
        T_list(inc_index+1,inc_index)=coh_tmm_bdata_list(nextstack_index,2);
    end
end
Ltilde=[1,-R_list(2,1);R_list(1,2),T_list(2,1)*T_list(1,2)-R_list(2,1)*R_list(1,2)]/T_list(1,2);
for i=2:num_inc_layers-1
    L=[1/P_list(i),0;0,P_list(i)]*[1,-R_list(i+1,i);R_list(i,i+1),T_list(i+1,i)*T_list(i,i+1)-R_list(i+1,i)*R_list(i,i+1)]/T_list(i,i+1);
    Ltilde=Ltilde*L;
end
T=1/Ltilde(1,1);
R=Ltilde(2,1)/Ltilde(1,1);
% we are only interested in calculating the R and T of incoherent case so
% the remaining part of tmm code is not included in this script.
end

