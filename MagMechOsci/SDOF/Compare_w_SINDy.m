clear all
close all
clc

%%
% Set true to use Ensemble SINDy
SINDy_ensemble=false; 

N_bootstrap=100;

%%
strain_idx=[5]; %
f_idx=12;
clock_beam_idx=4;
N_periods=200;



dim=length(strain_idx);
[freqs_up,x_POs_up,control_POs_up,t_POs_up]=get_POs_from_lvm('Forced_response_measurements/sweep_ampl_1p5_mag_0_up_0.lvm',strain_idx,f_idx,clock_beam_idx, N_periods,10^3);
freqs_up(end)=[];
x_POs_up(end)=[];
control_POs_up(end)=[];
t_POs_up(end)=[];

  x_POs_up(freqs_up<38 | freqs_up>41)=[];
  control_POs_up(freqs_up<38 | freqs_up>41)=[];
  t_POs_up(freqs_up<38 | freqs_up>41)=[];
  freqs_up(freqs_up<38 | freqs_up>41)=[];

  
 
 

[freqs_down,x_POs_down,control_POs_down,t_POs_down]=get_POs_from_lvm('Forced_response_measurements/sweep_ampl_1p5_mag_0_down_0.lvm',strain_idx,f_idx,clock_beam_idx,N_periods,10^3);

 freqs_down(end)=[];
x_POs_down(end)=[];
control_POs_down(end)=[];
t_POs_down(end)=[];

  x_POs_down(freqs_down<38 | freqs_down>41 )=[];
  control_POs_down(freqs_down<38 | freqs_down>41)=[];
  t_POs_down(freqs_down<38 | freqs_down>41)=[];
  freqs_down(freqs_down<38 | freqs_down>41)=[];
%x_POs_down=cellfun(@(c) [c(:,1) c(:,2) -c(:,3)],x_POs_down,'UniformOutput',false);
 
 




%%
x_POs=[x_POs_up; x_POs_down];
control_POs=[control_POs_up; control_POs_down];
t_POs=[t_POs_up; t_POs_down];
Om_vec=2*pi*[freqs_up  freqs_down];
%%
t_POs_shift=cell(size(t_POs));
for iter_POs=1:length(x_POs)
    control_hat=project_fmodes( t_POs{iter_POs},control_POs{iter_POs},1,Om_vec(iter_POs));
    phi0=angle(control_hat(3));
    t0=phi0./(Om_vec(iter_POs));
    t_POs_shift{iter_POs}=t_POs{iter_POs}+t0;
    
end

z_hat=NaN(length(Om_vec),dim,3);
for iter_Om=1:length(Om_vec)
    z_hat(iter_Om,:,:)=project_fmodes( t_POs_shift{iter_Om},x_POs{iter_Om},1,Om_vec(iter_Om));
end

figure 

plot(Om_vec./(2*pi),angle(z_hat(:,:,3)),'s')


%%
N_harm=10;
addpath(append(cd,'\sparsedynamics\sparsedynamics\utils'));
polyorder=3;
usesine=0;
n=2*dim;
f_order=2;
Theta_big=[];
dx=[];
for iter_Om=1:length(Om_vec)
    
     Om=Om_vec(iter_Om);
    x_hat=project_fmodes( t_POs{iter_Om},x_POs{iter_Om},N_harm,Om_vec(iter_Om));
    x_hat(N_harm+1)=0;
    
    v_hat=zeros(dim,2*N_harm+1);
    v_hat(1:N_harm)=x_hat(1:N_harm).*(1i*Om.*(-N_harm:1:-1));
    v_hat(N_harm+2:end)=x_hat(N_harm+2:end).*(1i.*Om.*(1:1:N_harm));
    
    a_hat=zeros(dim,2*N_harm+1);
    a_hat(1:N_harm)=x_hat(1:N_harm).*(-1).*(Om.*(-N_harm:1:-1)).^2;
    a_hat(N_harm+2:end)=x_hat(N_harm+2:end).*(-1).*(Om.*(1:1:N_harm)).^2;
    
    %x_POs_filt=real(x_hat*exp(1i.*Om*(-N_harm:N_harm).'*t_POs{iter_Om}.')).';
    %v_PO=real(v_hat*exp(1i.*Om*(-N_harm:N_harm).'*t_POs{iter_Om}.')).';
    %a_PO=real(a_hat*exp(1i.*Om*(-N_harm:N_harm).'*t_POs{iter_Om}.')).';
    T=2*pi/Om;
    t_proj=T*(0:0.001:10).';
    x_POs_filt=real(x_hat*exp(1i.*Om*(-N_harm:N_harm).'*t_proj.')).';
    v_PO=real(v_hat*exp(1i.*Om*(-N_harm:N_harm).'*t_proj.')).';
    a_PO=real(a_hat*exp(1i.*Om*(-N_harm:N_harm).'*t_proj.')).';
    
    x =[x_POs_filt v_PO];
    dx=[dx;v_PO a_PO];
    
    Theta = poolData(x,n,polyorder,usesine);
    Forcing=[];
    for iter_f_order=1:f_order
        %Forcing=[Forcing sin(iter_f_order*Om.*t_POs{iter_Om}) cos(iter_f_order*Om.*t_POs{iter_Om})];
        Forcing=[Forcing sin(iter_f_order*Om.*t_proj) cos(iter_f_order*Om.*t_proj)];

        %Forcing=[Forcing sin(iter_f_order*Om.*t_trans) cos(iter_f_order*Om.*t_trans)];
        %Forcing=[Forcing sin(iter_f_order*Om.*t_PO) cos(iter_f_order*Om.*t_PO)];
    end
    Theta_big =[Theta_big; Theta Forcing]; 
end
%%
lambda = 0.01;      % lambda is our sparsification knob.

if SINDy_ensemble==true
    idx_set=randi(length(Theta_big(:,1)),length(Theta_big(:,1)),N_bootstrap);
    %idx_set(:,1)=1:length(Theta(:,1));
    N_poly=nchoosek(polyorder+3-1,polyorder);
    Xi_ensemble=NaN(N_poly+2*f_order,n,N_bootstrap);
    
    for iter=1:N_bootstrap
        disp(['Bootstrap no:' num2str(iter)])
        Xi_ensemble(:,:,iter) = sparsifyDynamics(Theta_big(idx_set(:,iter),:),dx(idx_set(:,iter),:),lambda,n);
    end
    Xi=mean(Xi_ensemble,3);
else
    Xi = sparsifyDynamics(Theta_big,dx,lambda,n);
end
    

disp(['Eigenfrequency:  ' num2str(sqrt(-Xi(2,2))/(2*pi))])
disp(['Damping:  '        num2str( -Xi(3,2))])

disp(['q:  '          num2str( Xi(2,2))])
disp(['q_dot:  '      num2str( Xi(3,2))])

disp(['q^2:  '          num2str( Xi(4,2))])
disp(['q q_dot:  '      num2str( Xi(5,2))])
disp(['q_dot^2:  '      num2str( Xi(6,2))])


disp(['q^3:  '          num2str( Xi(7,2))])
disp(['q^2 q_dot:  '    num2str( Xi(8,2))])
disp(['q q_dot^2:  '    num2str( Xi(9,2))])
disp(['q_dot^3:  '      num2str( Xi(10,2))])

disp(['sin(Om*t):  '    num2str( Xi(11,2))])
disp(['cos(Om*t):  '    num2str( Xi(12,2))])
disp(['sin(2*Om*t):  '  num2str( Xi(13,2))])
disp(['cos(2*Om*t):  '  num2str( Xi(14,2))])

%%
z0=zeros(2*dim,1);
t_PO_sindy=cell(length(Om_vec),1);
x_PO_sindy=cell(length(Om_vec),1);
for iter_Om=1:length(Om_vec)
    Om=Om_vec(iter_Om);
    T=2*pi/(Om);
    
    RHS=@(t,x)[Xi(2,1)*x(1)+Xi(3,1)*x(2)+...
        Xi(4,1)*x(1)^2+Xi(5,1)*x(1)*x(2)+Xi(6,1)*x(2)^2+...
        Xi(7,1)*x(1)^3+Xi(8,1)*x(1)^2*x(2)+Xi(9,1)*x(1)*x(2)^2+Xi(10,1)*x(2)^3+...
        Xi(11,1)*sin(Om*t)+Xi(12,1)*cos(Om*t)+Xi(13,1)*sin(2*Om*t)+Xi(14,1)*cos(2*Om*t);
        Xi(2,2)*x(1)+Xi(3,2)*x(2)+...
        Xi(4,2)*x(1)^2+Xi(5,2)*x(1)*x(2)+Xi(6,2)*x(2)^2+...
        Xi(7,2)*x(1)^3+Xi(8,2)*x(1)^2*x(2)+Xi(9,2)*x(1)*x(2)^2+Xi(10,2)*x(2)^3+...
        Xi(11,2)*sin(Om*t)+Xi(12,2)*cos(Om*t)+Xi(13,2)*sin(2*Om*t)+Xi(14,2)*cos(2*Om*t)];
    
    % Transient simulations
    [t_trans, z_trans] = ode45(@(t,z)RHS(t,z), [0 500*T],z0); % Transients
    
%     figure
%     plot(t_trans,z_trans(:,1))
%     title(num2str(Om/(2*pi)))
    z0=z_trans(end,:)';
    % Steady State
    [t_PO, z_PO] = ode45(@(t,z)RHS(t,z), [0 10.*T], z0);
    
    if max(isnan(z0))==1
        z0=zeros(2*dim,1);
    end
    x_PO_sindy{iter_Om}=z_PO(:,1);
    t_PO_sindy{iter_Om}=t_PO;
    disp(['Progress: ' num2str(round(iter_Om/length(Om_vec)*100,2)) '%'])
end
%%

figure
[~,id]=max(Om_vec);
tmp_ampl=cell2mat(cellfun(@max,cellfun(@abs,x_POs,'UniformOutput',false),'UniformOutput',false));
plot(Om_vec(1:61)/(2*pi),tmp_ampl(1:61),'x')

hold on
plot(Om_vec(62:end)/(2*pi),tmp_ampl(62:end),'x')
Ampl_sindy=cell2mat(cellfun(@max,cellfun(@abs,x_PO_sindy,'UniformOutput',false),'UniformOutput',false));
 
plot(Om_vec(1:61)/(2*pi),Ampl_sindy(1:61),'x')

plot(Om_vec(62:end)/(2*pi),Ampl_sindy(62:end),'x')

legend('Experiment - Sweep Up','Experiment - Sweep Down', 'SinDy - Sweep Up','Sindy - Sweep Down')
xlabel('Excitation frequency in Hz')
ylabel('Amplitude |q^1| (V)')
grid on

%%
folder='decay_data';

cd(folder)
file_list=dir('*.lvm');
cd ..
strain_idx=3;

figure
for iter_file=1:1

data=lvm_import([folder '/' file_list(iter_file).name]);
t_start_idx=find(data.Segment1.data(:,11)==0,1)+1100;
t_end_idx=t_start_idx+3000;
t=data.Segment1.data(t_start_idx:t_end_idx,1);
t=t-t(1);
decay_sig=sgolayfilt(data.Segment1.data(t_start_idx:t_end_idx,strain_idx),10,41);
decay_sig=decay_sig-mean(decay_sig(t>2));


end

figure
plot(t,decay_sig)


%%

polyorder=5;
usesine=0;
n=2;


v_decay=gradient(decay_sig,t(2)-t(1));
a_decay=gradient(v_decay,t(2)-t(1));
% figure
% plot(t,decay_sig)
% hold on
% plot(t,v_decay)
% plot(t,a_decay)

%%
x=[decay_sig v_decay];
dx=[v_decay a_decay];
Theta = poolData(x,n,polyorder,usesine);
lambda = 0.01;      % lambda is our sparsification knob.


if SINDy_ensemble==true
    idx_set=randi(length(Theta(:,1)),length(Theta(:,1)),N_bootstrap);
    %idx_set(:,1)=1:length(Theta(:,1));
    N_poly=nchoosek(polyorder+3-1,polyorder);
    Xi_ensemble=NaN(N_poly,n,N_bootstrap);
    
    for iter=1:N_bootstrap
        
        Xi_ensemble(:,:,iter) = sparsifyDynamics(Theta(idx_set(:,iter),:),dx(idx_set(:,iter),:),lambda,n);
    end
    Xi=median(Xi_ensemble,3);
    
else
    Xi = sparsifyDynamics(Theta,dx,lambda,n);
end
poolDataLIST({'q','q_dot'},Xi,n,polyorder,usesine);

[t_sindy,x_sindy]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),t,x(1,:));  % approximate
figure
plot(t,decay_sig)
hold on
plot(t,x_sindy(:,1))
legend('Measurements','SinDy')
xlabel('time')
ylabel('Strain gauge signal (V)')

%%
f_ampl=0.131;

z0=zeros(2,1);


t_PO_sindy_decay=cell(length(Om_vec),1);
x_PO_sindy_decay=cell(length(Om_vec),1);
for iter_Om=1:length(Om_vec)
    Om=Om_vec(iter_Om);
    T=2*pi/(Om);
    RHS=@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine)+[0;f_ampl]*cos(Om*t);
    [t_trans, z_trans] = ode45(@(t,z)RHS(t,z), [0 500*T],z0); % Transients

%     figure
%     plot(t_trans,z_trans(:,1))
%     title(num2str(Om/(2*pi)))
%     z0=z_trans(end,:)';
    % Steady State
    [t_PO, z_PO] = ode45(@(t,z)RHS(t,z), [0 10.*T], z0);
    
    if max(isnan(z0))==1
        z0=zeros(2*dim,1);
        disp(['Reset initial conditions at frequency:' num2str(Om/(2*pi))])
    end
    x_PO_sindy_decay{iter_Om}=z_PO(:,1);
    t_PO_sindy_decay{iter_Om}=t_PO;
    disp(['Progress: ' num2str(round(iter_Om/length(Om_vec)*100,2)) '%'])
    
end
%%
figure
tmp_ampl=cell2mat(cellfun(@max,cellfun(@abs,x_POs,'UniformOutput',false),'UniformOutput',false));
plot(Om_vec/(2*pi),tmp_ampl/2,'x')
hold on
Ampl_sindy=cell2mat(cellfun(@max,cellfun(@abs,x_PO_sindy,'UniformOutput',false),'UniformOutput',false));
 
plot(Om_vec/(2*pi),Ampl_sindy/2,'x')

Ampl_sindy_decay=cell2mat(cellfun(@max,cellfun(@abs,x_PO_sindy_decay,'UniformOutput',false),'UniformOutput',false));
 
plot(Om_vec/(2*pi),Ampl_sindy_decay/2,'x')

legend('Experiment','SINDy with forced data', 'SINDy with decay data')
xlabel('Excitation frequency  (Hz)')
ylabel('Amplitude |q^1| (V)')
grid on

%%

function x_hat=project_fmodes(t_PO,x_PO,N_fmodes,Om)
x_hat=zeros(length(x_PO(1,:)),2*N_fmodes+1);
x_hat(:,N_fmodes+1)=mean(x_PO);

for iter_fmodes=1:N_fmodes
    %x_hat(:,N_fmodes+1+iter_fmodes)=1./(t_PO(end)-t_PO(1))*trapz(t_PO,x_PO.*exp(-1i*iter_fmodes*Om*t_PO)).';
    x_hat(:,N_fmodes+1+iter_fmodes)=(t_PO(2)-t_PO(1))./(t_PO(end)-t_PO(1))*sum(x_PO(2:end,:).*exp(-1i*iter_fmodes*Om*t_PO(2:end))).';
    x_hat(:,N_fmodes+1-iter_fmodes)=conj(x_hat(:,N_fmodes+1+iter_fmodes));
end

end


function [freqs, x_POs,control_POs,t_POs]=get_POs_from_lvm(src,strain_idx,f_idx,clock_beam_idx,N_periods,N_pts_per_T)

% data=lvm_import('nosignal_19042022.lvm');
% 
% base=mean(data.Segment1.data(:,strain_idx));
% 


data=lvm_import(src);

 
freq_idx=find(isnan(data.Segment1.data(:,f_idx))~=1);
freqs=data.Segment1.data(freq_idx,f_idx).';
freq_idx=freq_idx([1 1+find(diff(freqs)~=0)]);
freqs=freqs([1 1+find(diff(freqs)~=0)]);
x_POs=cell(length(freqs),1);
t_POs=cell(length(freqs),1);
control_POs=cell(length(freqs),1);

% mean(data.Segment1.data(1:freq_idx(2)-1,strain_idx))
for iter_freq=1:length(freq_idx)
    T=1/(data.Segment1.data(freq_idx(iter_freq),f_idx));
    
    Ts=T/N_pts_per_T;
    if iter_freq< length(freq_idx)
        t_end=freq_idx(iter_freq+1)-1;
    else
        t_end=length(data.Segment1.data(:,3)); 
    end
    t_start=find(data.Segment1.data(1:t_end,1)<data.Segment1.data(t_end,1)-N_periods*T,1,'last');
    
    
    
    control_beam_sig=data.Segment1.data(t_start:t_end,clock_beam_idx);
    
    control_beam_sig=control_beam_sig-mean(control_beam_sig);
    
    x_tmp=data.Segment1.data(t_start:t_end,strain_idx);
    
    x_tmp=x_tmp- mean(data.Segment1.data(1:freq_idx(2)-1,strain_idx));
    
    t_tmp=data.Segment1.data(t_start:t_end,1);
    
    t_interp=max(t_tmp)+(-N_periods*T:Ts:0);
    
    
    
  
    t_POs{iter_freq}=(t_interp-t_interp(1)).';
    x_POs{iter_freq}=reshape(interp1(t_tmp,x_tmp,t_interp,'spline','extrap'),length(t_interp),length(strain_idx));
    control_POs{iter_freq}=interp1(t_tmp,control_beam_sig,t_interp,'spline','extrap').';

    

end
end


