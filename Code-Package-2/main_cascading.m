%-------------------- load the S matrix elements for cascading ------------
num_modes_prop=init_data.num_modes_prop;
[N,M,ens_tot]=size(S21_array);
%------------------- Cascading rule ---------------------------------------
for ens_cas_count=1:no_of_thick_slabs
index=randperm(ens_tot);
%----------- first cascading ----------------------------------------------
S11_1=S11_array(:,:,index(1));
S21_1=S21_array(:,:,index(1));
S12_1=S12_array(:,:,index(1));
S22_1=S22_array(:,:,index(1));
%--------------------------------------------------------------------------
S11_2=S11_array(:,:,index(2));
S21_2=S21_array(:,:,index(2));
S12_2=S12_array(:,:,index(2));
S22_2=S22_array(:,:,index(2));
%--------------------------------------------------------------------------
S11_cas=S11_1 + S12_1*inv(eye(N,M)- S11_2*S22_1)*S11_2*S21_1;
S12_cas=S12_1*inv(eye(N,M)-S11_2*S22_1)*S12_2;
S21_cas=S21_2*inv(eye(N,M)-S22_1*S11_2)*S21_1;
S22_cas=S22_2 + S21_2*inv(eye(N,M)-S22_1*S11_2)*S22_1*S12_2;
 
%------------ From 3rd slab and beyond ------------------------------------   
   for cas_count=3:ens_tot
 
   S11_1=S11_cas;
   S12_1=S12_cas;
   S21_1=S21_cas;
   S22_1=S22_cas;

     S11_2=S11_array(:,:,index(cas_count));
     S21_2=S21_array(:,:,index(cas_count));
     S12_2=S12_array(:,:,index(cas_count));
     S22_2=S22_array(:,:,index(cas_count));
 %--------------------- cascading -----------------------------------------
   S11_cas=S11_1 + S12_1*inv(eye(N,M)- S11_2*S22_1)*S11_2*S21_1;
   S12_cas=S12_1*inv(eye(N,M)-S11_2*S22_1)*S12_2;
   S21_cas=S21_2*inv(eye(N,M)-S22_1*S11_2)*S21_1;
   S22_cas=S22_2 + S21_2*inv(eye(N,M)-S22_1*S11_2)*S22_1*S12_2;       
 %-------------------------------------------------------------------------
   S21_prop=S21_cas(1:num_modes_prop,1:num_modes_prop);  %% extract propagating TM
   S11_prop=S11_cas(1:num_modes_prop,1:num_modes_prop);  %% extract propagating RM
   S12_prop=S12_cas(1:num_modes_prop,1:num_modes_prop);  %% extract propagating TM
   S22_prop=S22_cas(1:num_modes_prop,1:num_modes_prop);  %% extract propagating RM

   %[cas_count trace(S21_prop'*S21_prop)/num_modes_prop]     
   end   
   
tau_avg=trace(S21_prop'*S21_prop)/num_modes_prop;      
% %---------------- Type 1 --------------------------------------------------
% x_hist=(0:0.01:1);
% fcount=zeros(1,length(x_hist));
% [U,Sigma,V] = svd(S21_prop); 
% tau=diag(Sigma.^2);
% hist_data=hist(tau,x_hist);
% area_f=trapz(x_hist,hist_data);
% fcount=(1./area_f).*hist_data;
% 
% 
% figure
% bar(x_hist,fcount);
% hold on
% T=linspace(sech(1/tau_avg)^2,0.99999,100000);
% ptrans= @(T) (tau_avg./(2.*T.*sqrt(1-T)));  %%% normalised eigen value density function 
% plot(T,ptrans(T),'LineWidth',2,'Color','r')
% %trapz(T,ptrans(T))  
% ylim([0 max(mean(fcount,1))+2])
% set(gca,'FontSize',18)
% xlabel('\tau')
% ylabel('p(\tau)')
% legend('Numerical modelling','R.M.T')
% %----------------- Type 2 : Normalised pdf --------------------------------
% tau_enh=tau./tau_avg;
% x_hist=(0:0.05:max(tau_enh));
% hist_data=hist(tau_enh,x_hist);
% area_f=trapz(x_hist,hist_data);
% fcount=(1./area_f).*hist_data;
% 
% figure
% bar(x_hist,fcount);
% set(gca,'FontSize',18)
% xlabel('$\frac{\tau}{<\tau>}$','Interpreter','Latex')
% ylabel('$p(\frac{\tau}{<\tau>})$','Interpreter','Latex')
%--------------------------------------------------------------------------
filename=sprintf('S21_cas_%d.mat',ens_cas_count);
save(filename,'S21_cas');
filename=sprintf('S12_cas_%d.mat',ens_cas_count);
save(filename,'S12_cas');
filename=sprintf('S11_cas_%d.mat',ens_cas_count);
save(filename,'S11_cas');
filename=sprintf('S22_cas_%d.mat',ens_cas_count);
save(filename,'S22_cas');
sprintf('Cascaded sample no %d with tau_avg = %.4f',ens_cas_count,tau_avg)
end
