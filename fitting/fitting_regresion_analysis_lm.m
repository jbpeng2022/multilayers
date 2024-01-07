function optimized_thickness=fitting_regresion_analysis_lm(d_fft,d_interval,nk,wavelength_all,theta,reflentance_experiment)
%% reproduce the train set
N=1000;%sampling point around the dfft
d_thickness_prediction_scope=linspace(d_fft-2*d_interval,d_fft+2*d_interval,N);
train_in=zeros(length(wavelength_all),N);
for i=1:N
   d=d_thickness_prediction_scope(i);
   R_T=tmm(d,nk,wavelength_all,theta);
   train_in(:,i)=R_T(:,1);
end
%% train the net
train_target=d_thickness_prediction_scope;
net=feedforwardnet(20, 'trainlm');
net.trainParam.showWindow=false;
net.trainParam.goal=1e-5;
net=train(net,train_in,train_target);
%% predict the output by using the net
target=reflentance_experiment;
optimized_thickness=net(target);
end