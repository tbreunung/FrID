clear all
close all
clc
% Download forced response masurements from https://github.com/mattiacenedese/BRBtesting
load('Forced_response_measurements/ForcedFrequencyResponses.mat')

plot(FRdata.ForceControl(1).frequency,FRdata.ForceControl(1).amplitude,'x')

hold on
plot(FRdata.ForceControl(2).frequency,FRdata.ForceControl(2).amplitude,'x')
plot(FRdata.ForceControl(3).frequency,FRdata.ForceControl(3).amplitude,'x')
plot(FRdata.ForceControl(4).frequency,FRdata.ForceControl(4).amplitude,'x')

force_nom=[0.2 0.5 1 1.5];

f0=1;

Om_mins=2*pi.*[79.2  79.2 79.2 79.2]/f0;
Om_maxs=2*pi.*[81.2  81.2 81.2 81.2]/f0;
ids=[1:4];
dim=1;

t_POs=cell(41,length(ids));
x_POs=cell(41,length(ids));
freqs=[];
ampls=[];
f_lvl=[];
for iter_ids=1:length(ids)
    freqs=[freqs FRdata.ForceControl(ids(iter_ids)).frequency.'./f0];
    ampls=[ampls FRdata.ForceControl(ids(iter_ids)).amplitude];
    phases=FRdata.ForceControl(ids(iter_ids)).phase;
    f_lvl=[f_lvl FRdata.ForceControl(ids(iter_ids)).force];
    for iter_POs=1:length(phases)
        T=1/freqs(iter_POs,iter_ids);
        dt=T/10^3;
        t_POs{iter_POs,iter_ids}=(0:dt:T).';
        x_POs{iter_POs,iter_ids}=ampls(iter_POs,iter_ids)/(2*pi*freqs(iter_POs,iter_ids))^2.*sin(2*pi*freqs(iter_POs,iter_ids).*t_POs{iter_POs,iter_ids}+phases(iter_POs));
    end
end
Om_vec=2*pi*freqs;
%%
x_POs_un={};
Om_vec_un=[];
t_POs_un={};
force_level=[];
id_f_lvl=[];
for iter_ids=1:length(ids)
    x_POs_tmp=x_POs(:,iter_ids);
    Om_tmp=Om_vec(:,iter_ids);
    t_POs_tmp=t_POs(:,iter_ids);
    x_POs_tmp=x_POs_tmp(:);
    Om_tmp=Om_tmp(:).';
    t_POs_tmp=t_POs_tmp(:);
    f_lvl_tmp=f_lvl(:,iter_ids);
    
    
    [Om_vec_un_tmp,x_POs_un_tmp,t_POs_un_tmp]=del_renun_orbits(Om_tmp,x_POs_tmp,t_POs_tmp,dim);
    
    tmp_ampl=cell2mat(cellfun(@max,cellfun(@abs,x_POs_un_tmp,'UniformOutput',false),'UniformOutput',false));
    
    tmp_ampl2=cell2mat(cellfun(@max,cellfun(@abs,x_POs_tmp,'UniformOutput',false),'UniformOutput',false));
    
    figure
    subplot(1,2,1)
    plot(Om_vec_un_tmp./(2*pi),sum(tmp_ampl,2),'x')
    title('Unique')
    subplot(1,2,2)
    plot(Om_tmp./(2*pi),sum(tmp_ampl2,2),'x')
    
    x_POs_un_tmp(Om_vec_un_tmp<Om_mins(iter_ids) | Om_vec_un_tmp>Om_maxs(iter_ids))=[];
    t_POs_un_tmp(Om_vec_un_tmp<Om_mins(iter_ids) | Om_vec_un_tmp>Om_maxs(iter_ids))=[];
    f_lvl_tmp(Om_vec_un_tmp<Om_mins(iter_ids) | Om_vec_un_tmp>Om_maxs(iter_ids))=[];
    Om_vec_un_tmp(Om_vec_un_tmp<Om_mins(iter_ids) | Om_vec_un_tmp>Om_maxs(iter_ids))=[];
    
    title('Repititions')
    
    x_POs_un=[x_POs_un;x_POs_un_tmp];
    t_POs_un=[t_POs_un;t_POs_un_tmp];
    Om_vec_un=[Om_vec_un Om_vec_un_tmp];
    force_level=[force_level; f_lvl_tmp ];
    id_f_lvl=[id_f_lvl; length(Om_vec_un)];
end

%%
N_harm=3;
x_POs_hat=zeros(length(Om_vec_un),1,2*N_harm+1);
for iter_Om=1:length(Om_vec_un)
    
    x_POs_hat(iter_Om,:,:)=project_fmodes( t_POs_un{iter_Om},x_POs_un{iter_Om},N_harm,Om_vec_un(iter_Om));
    
end





%%


figure
tl = tiledlayout(2,2);
nexttile(tl)
tmp_ampl=cell2mat(cellfun(@max,cellfun(@abs,x_POs_un,'UniformOutput',false),'UniformOutput',false));
%
plot(Om_vec_un./(2*pi),tmp_ampl,'x')

title('max displacement')


for iter_harm=1:N_harm
    nexttile(tl)
    plot(Om_vec_un./(2*pi),abs(x_POs_hat(:,:,N_harm+1+iter_harm)),'x')
    
    
    title([ num2str(iter_harm) 'th harmonic'])
    xlabel('Frequency in Hz')
    ylabel('Amplitude  in V')
    
end
title(tl,'Response ')



%%


norm_c=(max(cell2mat(cellfun(@max,cellfun(@abs,x_POs_un,'UniformOutput',false),'UniformOutput',false))));

x_POs_un=gdivide(x_POs_un,norm_c);



%%





z_hat=NaN(length(Om_vec_un),dim,3);
for iter_Om=1:length(Om_vec_un)
    z_hat(iter_Om,:,:)=project_fmodes( t_POs_un{iter_Om},x_POs_un{iter_Om},1,Om_vec_un(iter_Om));
end

figure

plot(Om_vec_un./(2*pi),angle(z_hat(:,:,3)),'s')
%%

noise_lvl=noise_estimator(x_POs_un,t_POs_un,Om_vec_un,500,5000*2*pi);

%%
poly_order=5;
N_fmodes=1;
f_order=1;




[Abig, Bbig]=Build_big_mat(x_POs_un,t_POs_un,force_level,N_fmodes,Om_vec_un,dim,poly_order,f_order,noise_lvl);
%%
ws=zeros(length(Abig(:,1)),1);

figure

for iter_flvl=1:length(id_f_lvl)
    if iter_flvl==1
        Om_tmp=repmat(Om_vec_un(1:id_f_lvl(iter_flvl)),2*N_fmodes,1);
        
        ws(1:2*N_fmodes*id_f_lvl(iter_flvl))=1/max(abs(Abig(1:2*N_fmodes*id_f_lvl(iter_flvl),1)));
        ampl_tmp=ws(1:2*N_fmodes*id_f_lvl(iter_flvl)).*abs(Abig(1:2*N_fmodes*id_f_lvl(iter_flvl),1));
        ws(1:2*N_fmodes*id_f_lvl(iter_flvl))=ws(1:2*N_fmodes*id_f_lvl(iter_flvl))./id_f_lvl(iter_flvl); %
    else
        Om_tmp=repmat(Om_vec_un(id_f_lvl(iter_flvl-1)+1:id_f_lvl(iter_flvl)),2*N_fmodes,1);
        ws(2*N_fmodes*id_f_lvl(iter_flvl-1)+1:2*N_fmodes*id_f_lvl(iter_flvl))=1/max(abs(Abig(2*N_fmodes*id_f_lvl(iter_flvl-1)+1:2*N_fmodes*id_f_lvl(iter_flvl),1)));
        
        ampl_tmp=ws(2*N_fmodes*id_f_lvl(iter_flvl-1)+1:2*N_fmodes*id_f_lvl(iter_flvl)).*abs(Abig(2*N_fmodes*id_f_lvl(iter_flvl-1)+1:2*N_fmodes*id_f_lvl(iter_flvl),1));
        ws(2*N_fmodes*id_f_lvl(iter_flvl-1)+1:2*N_fmodes*id_f_lvl(iter_flvl))=ws(2*N_fmodes*id_f_lvl(iter_flvl-1)+1:2*N_fmodes*id_f_lvl(iter_flvl))./(id_f_lvl(iter_flvl)-id_f_lvl(iter_flvl-1));%
    end
    plot(Om_tmp(:),ampl_tmp,'x')
    hold on
end
ws=ws./min(ws);



%%
tols=0.0015;
poly_str=write_polys(dim,poly_order,f_order);
pars_exp=zeros(size(poly_str));
err_rel=zeros(dim,1);
poly_str_tmp=poly_str(:,1);
poly_str_tmp=cellfun(@(c) c(8:end) ,poly_str_tmp,'UniformOutput',false);
for iter_dim=1:dim
    [pars_exp(:,iter_dim), ~,err_rel(iter_dim)]=FourierID(ws.*Abig,ws.*Bbig(iter_dim,:).',noise_lvl,tols(iter_dim),poly_str_tmp); %([1:83 250 251])([1:83 250 251])
end

%%
for iter_dim=1:dim
    disp('---------------------')
    disp('ID result:')
    disp(['Relative error:   ' num2str(err_rel(iter_dim))] )
    disp('---------------------')
    for iter_pars=1:length(pars_exp)
        if ~isnan(pars_exp(iter_pars,iter_dim))
            disp(['Coefficient: ' poly_str{iter_pars,iter_dim}  '     Value: ' num2str(pars_exp(iter_pars,iter_dim))])
        end
    end
    disp('---------------------')
    disp('Deleted coefficients')
    for iter_pars=1:length(pars_exp)
        if isnan(pars_exp(iter_pars))
            disp(['Coefficient: ' poly_str{iter_pars,iter_dim}  ])
        end
    end
    disp('---------------------')
    
end



%%
ampl_id= pars_exp(end-2*f_order+1:end,:);
pars_sim=real(pars_exp);
pars_sim(end-2*f_order+1:end,:)=[];
pars_sim=pars_sim.';

ampl_id=reshape(ampl_id,2*f_order,dim).';

ampl_id(:,end-f_order+2:end+1)=ampl_id(:,end-f_order+1:end);
ampl_id(:,f_order+1)=0;



pars_sim(isnan(pars_sim))=0;
ampl_id(isnan(ampl_id))=0;


IC=zeros(1,2*dim);
Om_step=(max(Om_vec_un)-min(Om_vec_un))/100;
Om_sim=(min(Om_vec_un):Om_step:max(Om_vec_un));
T_sim=200*2*pi/mean(Om_sim);

f_levels=force_nom;
x_hat_id=zeros(length(f_levels),length(Om_sim),dim,2*N_harm+1);
parfor iter_f_lvl=1:length(f_levels)
    [~,~, ~,x_tmp]=my_sweep_sim(Om_sim,pars_sim,poly_order,dim,f_levels(iter_f_lvl).*ampl_id,T_sim,T_sim,0.001,N_harm,IC);
    x_hat_id(iter_f_lvl,:,:,:)=x_tmp;
end
%%
z_hat=NaN(length(Om_vec_un),dim,2*N_harm+1);
for iter_Om=1:length(Om_vec_un)
    z_hat(iter_Om,:,:)=project_fmodes( t_POs_un{iter_Om},x_POs_un{iter_Om},N_harm,Om_vec_un(iter_Om));
    
end

%%


%%

figure
iter_harm=1;


for iter_f_lvl=1:length(f_levels)
    legend('AutoUpdate','Off')
    plot(freqs(:,iter_f_lvl),ampls(:,iter_f_lvl)/2,'x')
    hold on
    legend('AutoUpdate','On')
    set(gca,'ColorOrderIndex',iter_f_lvl)
    [x_hat_plt,Om_plot]=id_jump(reshape(x_hat_id(iter_f_lvl,:,:,:),[length(Om_sim) dim 2*N_harm+1]),Om_sim,N_harm+2);
    plot(Om_plot./(2*pi),(Om_plot.').^2.*sum(squeeze(abs(x_hat_plt(:,:,N_harm+1+iter_harm))),2).*norm_c,'Linewidth',2,'Displayname',['f=' num2str(f_levels(iter_f_lvl))])%,'Color',[0.8 0.3 0.3]
end
grid on
xlabel('Frequency in Hz')
ylabel('Acceleration amplitude in m/s^2')


%%
function [Abig, Bbig]=Build_big_mat(z_POs,t_POs,force_level,N_fmodes,Om_vec,dim,poly_order,f_order,noise_lvl)


disp('-----------------------------')
disp('Builing matrix for Regression')
disp('-----------------------------')
poly_str=write_polys(dim,poly_order,f_order);
N_orbits=length(z_POs);

Abig3D=zeros(2*N_fmodes,length(poly_str),N_orbits);
Bbig3D=zeros(dim,2*N_fmodes,N_orbits);
for iter_POs=1:N_orbits
    [A0,B0]=build_mat(z_POs{iter_POs},t_POs{iter_POs},force_level(iter_POs),N_fmodes,Om_vec(iter_POs),dim,poly_order,f_order,noise_lvl);
    % delte mean rows
    A0(N_fmodes+1,:)=[];
    % delete zeroth forcing coefficient
    A0(:,end-f_order)=[];
    % delete zeroth acceleration coefficient
    B0(:,N_fmodes+1)=[];
    Abig3D(:,:,iter_POs)=A0;
    Bbig3D(:,:,iter_POs)=B0;
    
    disp(['Progress: ' num2str(round(iter_POs/N_orbits*100,2)) '%'])
end
Abig=reshape(permute(Abig3D,[1 3 2]),2*N_fmodes*N_orbits,length(poly_str));
Bbig=reshape(Bbig3D,dim,2*N_fmodes*N_orbits);


disp('-----------------------------')
disp('Matrix Build')
disp('-----------------------------')
end

function [pars, res,rel_err]=FourierID(Abig,Bbig,noise_lvl,tol,poly_str)
disp('-----------------------------')
disp('Start fitting')
disp('-----------------------------')

del_idx=find(vecnorm(Abig)<noise_lvl);

accept_idx=min(setdiff(1:length(Abig(1,:)),del_idx));
del_idx=[];
for iter_pars=accept_idx+1:length(Abig(1,:))
    
    P=Abig(:,accept_idx)*inv(Abig(:,accept_idx)'*Abig(:,accept_idx))*Abig(:,accept_idx)';
    if  norm(Abig(:,iter_pars))<noise_lvl  || norm(P*Abig(:,iter_pars))/norm(Abig(:,iter_pars))>0.99
        
        del_idx=[del_idx iter_pars];
        
    else
        
        accept_idx=[accept_idx iter_pars];
        
    end
end
%


if ~isempty(del_idx)
    disp('Renundant Data detected!')
    disp(['Indices: ' num2str(del_idx)  ])
    
    disp(strjoin(['Polynomials: ' poly_str(del_idx).'],'   '))
end


accept_idx=setdiff(1:length(Abig(1,:)),del_idx);


A_reduce=Abig;
A_reduce(:,del_idx)=[];
pars=-A_reduce\Bbig;

res=norm(A_reduce*pars+Bbig);

if res/norm(Bbig)<tol
    disp(['Relative errror:' num2str(res/norm(Bbig)) ' less than ' num2str(tol*100) '%'])
    disp('Sparsify')
    
    chck=false;
    final_accept_idx=[];
    while chck==false
        
        res_tmp=zeros(length(accept_idx),1);
        pars_tmp=zeros(length(final_accept_idx)+1,length(accept_idx));
        
        for iter_idx=1:length(accept_idx)
            A_reduce=Abig(:,[final_accept_idx  accept_idx(iter_idx)]);
            
            pars_tmp(:,iter_idx)=-A_reduce\Bbig;
            
            res_tmp(iter_idx)=norm(A_reduce*pars_tmp(:,iter_idx)+Bbig);
            
        end
        [res_min, min_id]=min(res_tmp);
        if res_min/norm(Bbig)>tol
            disp(['Identified index: ' num2str(accept_idx(min_id))])
            disp(['Identified Polynomial:   ' poly_str{accept_idx(min_id)}])
            disp(['Relative error: ' num2str(res_min/norm(Bbig))])
            final_accept_idx=[final_accept_idx accept_idx(min_id)];
            accept_idx(min_id)=[];
            pars=pars_tmp(:,min_id);
            res=res_min;
        else
            disp(['Identified index: ' num2str(accept_idx(min_id))])
            disp(['Identified Polynomial:   ' poly_str{accept_idx(min_id)}])
            disp(['Relative error: ' num2str(res_min/norm(Bbig))])
            final_accept_idx=[final_accept_idx accept_idx(min_id)];
            accept_idx(min_id)=[];
            pars=pars_tmp(:,min_id);
            res=res_min;
            disp(['Sparse model with ' num2str(tol*100) '% relative error fitted'])
            chck=true;
        end
        if isempty(accept_idx)
            disp(['Model cannot be fitted with reaching a relative error of ' num2str(tol*100) '%'])
            chck=true;
        end
        
        
    end
else
    disp(['Relative errror:' num2str(res/norm(Bbig)) ' larger than specified (' num2str(tol*100) '%).'])
    disp(['No spare model fitted. Consider increasing tolercane or generating cleaner input signals.'])
    final_accept_idx=accept_idx;
end

del_idx=setdiff(1:length(Abig(1,:)),final_accept_idx);

[~,ids]=sort([final_accept_idx del_idx].');

pars=[  pars; NaN(length(del_idx),1)];

pars=pars(ids);
rel_err=res./norm(Bbig);
end

function str=write_polys(dim,poly_order,f_order)

cords=cell(2*dim,1);
for iter_dim=1:dim
    cords{iter_dim}=['x_' num2str(iter_dim)];
end
for iter_dim=1:dim
    cords{dim+iter_dim}=['v_' num2str(iter_dim)];
end

str=cords;



for iter_poly=2:poly_order
    str=[str; join(nmultichoosek(cords, iter_poly))];
end
num_poly=length(str);
str=repmat(str,1,dim);
prefix={};

for iter_dim=1:dim
    prefix=[prefix  repmat({['eq ' num2str(iter_dim) ':   ']},num_poly,1)];
end
str=strcat(prefix,str);

f_str=cell(2*f_order,1);
for iter_f=1:f_order
    f_str{f_order+iter_f}=['f^' num2str(iter_f) '_+'];
    f_str{f_order-iter_f+1}=['f^' num2str(iter_f) '_-'];
end
f_str=repmat(f_str,1,dim);

prefix=[];
for iter_dim=1:dim
    prefix=[prefix repmat({['eq ' num2str(iter_dim) ':   ']},2*f_order,1)];
end

f_str=strcat(prefix,f_str);

str=[str; f_str];
end


function [A,B]=build_mat(x_PO,t_PO,f_lvl,N_fmodes,Om,dim,poly_order,f_order,noise_floor)

if noise_floor>0
    x_hat=project_fmodes(t_PO,x_PO,N_fmodes,Om);
    x_hat(abs(x_hat)<noise_floor)=0;
    v_hat=(1i*Om.*(-N_fmodes:1:N_fmodes)).*x_hat;
    v_hat(abs(v_hat)<noise_floor)=0;
    x_PO=(x_hat*exp(1i.*Om*(-N_fmodes:N_fmodes).'*t_PO.')).';
    v_PO=(v_hat*exp(1i.*Om*(-N_fmodes:N_fmodes).'*t_PO.')).';
end
z_PO=[x_PO v_PO];

num_polys=0;
for iter_poly=1:poly_order
    num_polys=num_polys+nchoosek(2*dim+iter_poly-1,iter_poly);
end
polys_hat=zeros(2*N_fmodes+1,num_polys);
ids=1;
for iter_poly=1:poly_order
    idx=nmultichoosek(1:2*dim, iter_poly);
    for iter_idx=1:length(idx(:,1))
        sig=prod(z_PO(:,idx(iter_idx,:)),2);
        polys_hat(:,ids)=project_fmodes(t_PO,sig,N_fmodes,Om);
        ids=ids+1;
    end
end
polys_hat(abs(polys_hat)<noise_floor)=0;

A=polys_hat;
B=-Om^2*(-N_fmodes:N_fmodes).^2.*x_hat;
I_f=f_lvl.*[zeros((N_fmodes-f_order),(2*f_order+1));eye((2*f_order+1));zeros((N_fmodes-f_order),(2*f_order+1))];


A=[A I_f];

end

function [z_POs,t_POs, z_ampl,x_hat]=my_sweep_sim(Om_vec,pars,poly_order,dim,ampls,t_wait,N_periods_ss,t_sample_ss,N_harm,IC)

opts = odeset();

z_POs=cell(length(Om_vec),1);
t_POs=cell(length(Om_vec),1);
z_ampl=NaN(length(Om_vec),dim);
x_hat=NaN(length(Om_vec),dim,N_harm*2+1);

t0=0;
for iter_Om=1:length(Om_vec)
    Om_cur=Om_vec(iter_Om);
    T_cur=2*pi/Om_cur;
    
    [t_tmp, z] = ode45(@(t,z)NL_sys(t,z,pars,poly_order,dim,ampls,Om_cur),t0+[0 t_wait], IC ,opts); %;zeros(4*dim^2,1)
    
    t0= rem(t_tmp(end),2*pi/Om_cur);
    
    [t_tmp, z] = ode45(@(t,z)NL_sys(t,z,pars,poly_order,dim,ampls,Om_cur),t0+T_cur*[0:t_sample_ss:N_periods_ss], z(end,:) ,opts);
    
    t0= rem(t_tmp(end),2*pi/Om_cur);
    t_POs{iter_Om}=t_tmp;
    z_POs{iter_Om}=z(:,1:dim);
    
    IC=z(end,:);
    z_ampl(iter_Om,:)=max(abs(z(:,1:dim)));
    x_hat(iter_Om,:,:)=project_fmodes(t_tmp,z(:,1:dim),N_harm,Om_cur);
    
    
    disp(['Progress: ' num2str(round(iter_Om/length(Om_vec)*100,2)) '%'])
end
end

function noise_lvl=noise_estimator(z_POs,t_POs,Om_vec,N_smpls,Om_max)
idx_selected=randi(length(z_POs),N_smpls,1);
noise_ampls=zeros(length(z_POs{1}(1,:)),2*N_smpls);


for iter_smpls=1:N_smpls
    
    Om_found=false;
    while Om_found==false
        Oms_selected=Om_max*rand(1);
        if mod(Oms_selected,Om_vec(idx_selected(iter_smpls)))>0.1*Om_vec(idx_selected(iter_smpls))
            Om_found=true;
        end
    end
    
    x_hat_tmp=project_fmodes(t_POs{idx_selected(iter_smpls)},z_POs{idx_selected(iter_smpls)},1,Oms_selected);
    noise_ampls(:,2*iter_smpls-1:2*iter_smpls)=abs(x_hat_tmp(:,1:2:3));
end
noise_lvl=2*mean(noise_ampls,'all');
end


function x_hat=project_fmodes(t_PO,x_PO,N_fmodes,Om)
x_hat=zeros(length(x_PO(1,:)),2*N_fmodes+1);
x_hat(:,N_fmodes+1)=mean(x_PO);

for iter_fmodes=1:N_fmodes
    x_hat(:,N_fmodes+1+iter_fmodes)=(t_PO(2)-t_PO(1))./(t_PO(end)-t_PO(1))*sum(x_PO(2:end,:).*exp(-1i*iter_fmodes*Om*t_PO(2:end))).';
    x_hat(:,N_fmodes+1-iter_fmodes)=conj(x_hat(:,N_fmodes+1+iter_fmodes));
end

end

function [x_hat_plt,Om_plot]=id_jump(x_hat_with_jump,Om_with_jump,idx_for_jump_detection)

x_hat_plt=x_hat_with_jump;
Om_plot=Om_with_jump;

dim2=length(x_hat_with_jump(1,:,1));
dim3=length(x_hat_with_jump(1,1,:));
jump_idx=find(abs(diff(sum(abs(x_hat_with_jump(:,:,idx_for_jump_detection)),2)))>0.1);

jump_idx=jump_idx+[0:(length(jump_idx)-1)].';
for iter_jump=1:length(jump_idx)
    x_hat_plt=[x_hat_plt(1:jump_idx(iter_jump),:,:); NaN(1,dim2,dim3);x_hat_plt(jump_idx(iter_jump)+1:end,:,:) ];
    Om_plot=[Om_plot(1:jump_idx(iter_jump)) NaN Om_plot(jump_idx(iter_jump)+1:end)];
    
end
end



function x_p = NL_sys(t,x,pars,poly_order,dim,ampls,Om)
tol=10^-6;
f=sum(ampls.*exp(1i*Om.*(( 1:length(ampls))-floor((length(ampls)+1)/2))*t),2);
if imag(f)>tol
    disp('img forcing')
    
end
f=real(f);

RHS=zeros(dim,1);


polys=[x.'];
for iter_poly=2:poly_order
    polys=[polys  prod(nmultichoosek(x,iter_poly).')];
    
end
for iter_dim=1:dim
    RHS(iter_dim)=dot(pars(iter_dim,:),polys);
end
x_p=zeros(2*dim,1);
x_p(1:dim,1) = x(dim+1:2*dim);
x_p(dim+1:2*dim,1) = -RHS+f;
end


function y=BeamID_system(t,x,p)

x1 = x(1,:);
x2 = x(2,:);
p1 = p(1,:);
p2 = p(2,:);


d = p(3,:);

w0 = sqrt(56979.5144);
k3 = 15572.8828;
k2 =  10148.5048;

y(1,:) = x2;
y(2,:) = -d.*x2-w0.^2.*x1-k3.*x1.^3-k2.*x1.^2+p2.*cos(2*pi./p1.*t); %
end

function combs = nmultichoosek(values, k)
%// Return number of multisubsets or actual multisubsets.
% thanks to https://stackoverflow.com/questions/28284671/generating-all-combinations-with-repetition-using-matlab
if numel(values)==1
    n = values;
    combs = nchoosek(n+k-1,k);
else
    n = numel(values);
    combs = bsxfun(@minus, nchoosek(1:n+k-1,k), 0:k-1);
    combs = reshape(values(combs),[],k);
end
end

function [Oms,z_PO,t_PO]=del_renun_orbits(Oms,z_PO,t_PO,dim)
freqs_unique=unique(Oms);
freq_del=[];
for iter_freqs=1:length(freqs_unique)
    freq_idx=find(freqs_unique(iter_freqs)==Oms);
    x_ampl=zeros(length(freq_idx),dim);
    for iter_id=1:length(freq_idx)
        x_hat=project_fmodes(t_PO{freq_idx(iter_id)},z_PO{freq_idx(iter_id)},1,Oms(freq_idx(iter_id)));
        x_ampl(iter_id,:)=abs(x_hat(:,3)).';
    end
    del_id=[];
    for iter_ampl=1:length(x_ampl(:,1))-1
        ampl_id=find(norm(x_ampl(iter_ampl,:)-x_ampl(iter_ampl+1:end,:))<0.2*norm(x_ampl(iter_ampl,:)));
        del_id=unique([del_id iter_ampl+ampl_id]);
    end
    freq_del=[freq_del freq_idx(del_id)];
end
z_PO(freq_del)=[];

t_PO(freq_del)=[];%;
Oms(freq_del)=[];
end

function [freqs, x_POs,control_POs,t_POs]=get_POs_from_lvm(src,strain_idx,f_idx,clock_beam_idx,N_periods,N_pts_per_T)


data=lvm_import(src);


freq_idx=find(isnan(data.Segment1.data(:,f_idx))~=1);
freqs=data.Segment1.data(freq_idx,f_idx).';
freq_idx=freq_idx([1 1+find(diff(freqs)~=0)]);
freqs=freqs([1 1+find(diff(freqs)~=0)]);
x_POs=cell(length(freqs),1);
t_POs=cell(length(freqs),1);
control_POs=cell(length(freqs),1);

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