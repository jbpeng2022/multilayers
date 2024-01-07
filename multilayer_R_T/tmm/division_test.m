clc
clear
n_list=[1,2.2+1i*0.5,2.2,2.3,2.6,5.5+1i*0.7,1];
d_list=[inf,150,140,1000,1500,160,inf];
c_list={'i','c','c','i','i','c','i'};

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

% num_stacks
% num_inc_layers
% inc_from_all
% stack_from_all
% inc_index
% stack_index
% all_from_stack
% all_from_inc
% stack_d_list
% stack_n_list
% stack_from_inc

% false
% the format of c_list is cell rather than array so the way of calling the
% element is {}.
% 
% 
% 
% 