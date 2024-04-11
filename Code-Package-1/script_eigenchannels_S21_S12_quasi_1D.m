%------------------------- Eigenchannel decomposition ---------------------
% plot options are 'real_part','imaginary_part','magnitude'
%plot_type='real_part';  
%plot_type='imaginary_part';  
plot_type='magnitude';  
%------------------ Part 1 : Wave incidence from left to right ------------
[U_S21,Sigma_S21,V_S21] = svd(S21_pr_pr);% Eigenchannel modes estimated only for propagating modes 
tau_S21=diag(Sigma_S21).^2;             % Transmission of eigenchannels
tau_avg_S21=mean(tau_S21);               % Average transmission
%--- Estimate probability density of eigenchannel transmission coefficients ----------------------
no_of_bins=30;
[tau_hist_S21,x_hist] = histcounts(tau_S21,linspace(0,1,no_of_bins),'Normalization', 'pdf');
mid_point = x_hist(1:end-1)+mean(diff(x_hist))/2;   
dtau=(mid_point(2)-mid_point(1));
[sum(tau_hist_S21.*dtau) 1]  % Integrating the pdf
[tau_avg_S21 trace(S21_pr_pr'*S21_pr_pr)./init_data.num_modes_prop] % Average transmission

figure('position',[100 100 1500 600])
subplot(1,2,1)
plot(1:init_data.num_modes_prop,tau_S21,'-*')
xlabel('$Eigenchannel~no$','Interpreter','Latex')
ylabel('$\tau$','Interpreter','Latex')
title('$Transmission~eigenchannels~of~S_{21}$','Interpreter','Latex')
set(gca,'FontSize',18);
subplot(1,2,2)
bar(mid_point,tau_hist_S21);
xlabel('$\tau$','Interpreter','Latex')
ylabel('$p(\tau)$','Interpreter','Latex')
title({'$P.D.F~of~eigenchannel~transmission$','$~coefficient~\tau~of~S_{21}$'},'Interpreter','Latex')
set(gca,'FontSize',18) 
if tau_avg_S21 < 0.25  % If in the diffusion regime, use the bimodal pdf 
    % from the random matrix theory for comparison
    tau_val=linspace(sech(1/tau_avg_S21)^2,0.99999,100000);
    ptrans= @(tau_val,tau_avg) (tau_avg./(2.*tau_val.*sqrt(1-tau_val))); 
    % Bimodal density function 
    hold on
plot(tau_val,ptrans(tau_val,tau_avg_S21),'LineWidth',2,'Color','r')
trapz(tau_val,ptrans(tau_val,tau_avg_S21))  
set(gca, 'YScale', 'log')
legend('$From~modelling$','$From~R.M.T$','Interpreter','Latex')
end
set(gca,'FontSize',18);
side_description={sprintf('$\\tau_{max}=%.4f$',tau_S21(1)), ...
    sprintf('$\\tau_{min}=%.3e$',tau_S21(end)),...
    sprintf('$<\\tau>=%.4f$',mean(tau_S21))};
annotation('textbox', [0.005, 0.9, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on');   
%--------------------------------------------------------------------------
c_inc=V_S21(:,1);           % Maximally transmitting (First) singular vector
c_trans=S21_pr_pr*c_inc;     % Transmitted eigenchannel wavefront
sprintf('Transmission = %f',sum(abs(c_trans(1:init_data.num_modes_prop).^2)))
% S21_pr_pr contains only propagating modes, hence a truncated transmission 
% matrix containing both propagating and evanescent waves on the outgoing side
% and only propagating modes on the incident side is defined for
% visualisation purposes for applying the boundary integral equation.
S21_for_eigen=S21(:,1:init_data.num_modes_prop);
S11_for_eigen=S11(:,1:init_data.num_modes_prop);
% maximally transmitting eigenchannel
field_visualisation_and_comparison_LR(c_inc,Gij_LR,S21_for_eigen,S11_for_eigen,init_data,plot_type)
side_description={'$First (Maximally $','$transmitting)$','$eigenchannel~of~S_{21}$'};
annotation('textbox', [0.005, 0.95, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on');   
%-------- An intermediate eigenchannel ------------------------------------
eig_channel_mode_no=floor(2*(init_data.num_modes_prop)/3)+1;
c_inc=V_S21(:,eig_channel_mode_no); % transmitting singular vector
c_trans=S21_pr_pr*c_inc; % Transmitted eigenchannel wavefront
%c_refl=S11_pr_pr*c_inc; % Transmitted eigenchannel wavefront
sprintf('Transmission = %f',sum(abs(c_trans(1:init_data.num_modes_prop).^2)))
field_visualisation_and_comparison_LR(c_inc,Gij_LR,S21_for_eigen,S11_for_eigen,init_data,plot_type)
side_description={sprintf('$~%d ^{th}~eigenchannel$',eig_channel_mode_no),'$of~S_{21}$',};
annotation('textbox', [0.005, 0.95, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on'); 

%---------- Minimally transmitting eigenchannel ---------------------------
c_inc=V_S21(:,end);     % Minimally transmitting singular vector
c_trans=S21_pr_pr*c_inc; % Transmitted eigenchannel wavefront
c_refl=S11_pr_pr*c_inc; % Transmitted eigenchannel wavefront
sprintf('Transmission = %f',sum(abs(c_trans(1:init_data.num_modes_prop).^2)))
field_visualisation_and_comparison_LR(c_inc,Gij_LR,S21_for_eigen,S11_for_eigen,init_data,plot_type)
side_description={sprintf('$%d^{th} (Minimally $',size(V_S21,2)),'$transmitting)$','$eigenchannel~of~S_{21}$'};
annotation('textbox', [0.005, 0.95, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on'); 


%------------------ Part 2 : Wave incidence from right to left ------------
%---------- Maximally transmitting eigenchannel ---------------------------
[U_S12,Sigma_S12,V_S12] = svd(S12_pr_pr);% Eigenchannel modes estimated only for propagating modes 
tau_S12=diag(Sigma_S12).^2;    % Transmission of eigenchannels
c_inc=V_S12(:,1);              % Maximally transmitting (First) singular vector
c_trans=S12_pr_pr*c_inc;     % Transmitted eigenchannel wavefront
sprintf('Transmission = %f',sum(abs(c_trans(1:init_data.num_modes_prop).^2)))
S12_for_eigen=S12(:,1:init_data.num_modes_prop);
S22_for_eigen=S22(:,1:init_data.num_modes_prop);
field_visualisation_and_comparison_RL(c_inc,Gij_RL,S12_for_eigen,S22_for_eigen,init_data,plot_type)
side_description={'$First (Maximally $','$transmitting)$','$eigenchannel~of~S_{12}$'};
annotation('textbox', [0.005, 0.95, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on');   

%-------- An intermediate eigenchannel ------------------------------------
eig_channel_mode_no=floor(2*(init_data.num_modes_prop)/3)+1;
c_inc=V_S12(:,eig_channel_mode_no);         % Minimally transmitting singular vector
c_trans=S12_pr_pr*c_inc;     % Transmitted eigenchannel wavefront
sprintf('Transmission = %f',sum(abs(c_trans(1:init_data.num_modes_prop).^2)))
S12_for_eigen=S12(:,1:init_data.num_modes_prop);
S22_for_eigen=S22(:,1:init_data.num_modes_prop);
field_visualisation_and_comparison_RL(c_inc,Gij_RL,S12_for_eigen,S22_for_eigen,init_data,plot_type)
side_description={sprintf('$~%d ^{th}~eigenchannel$',eig_channel_mode_no),'$of~S_{12}$',};
annotation('textbox', [0.005, 0.95, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on'); 
%---------- Minimally transmitting eigenchannel ---------------------------
c_inc=V_S12(:,end);         % Minimally transmitting singular vector
c_trans=S12_pr_pr*c_inc;     % Transmitted eigenchannel wavefront
sprintf('Transmission = %f',sum(abs(c_trans(1:init_data.num_modes_prop).^2)))
S12_for_eigen=S12(:,1:init_data.num_modes_prop);
S22_for_eigen=S22(:,1:init_data.num_modes_prop);
field_visualisation_and_comparison_RL(c_inc,Gij_RL,S12_for_eigen,S22_for_eigen,init_data,plot_type)
side_description={sprintf('$%d^{th} (Minimally $',size(V_S21,2)),'$transmitting)$','$eigenchannel~of~S_{12}$'};
annotation('textbox', [0.005, 0.95, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on'); 
