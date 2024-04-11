%-------------------- 2D transport problem --------------------------------
% Estimate S21_lead for creating a finite sized incident modes with
%      transverse size $D_{in}< W$
% Analysis can be reduced to a diffraction problem where the light travels
%      from a smaller waveguide to a larger one. 
% As there is no flux loss, S21_lead has to be a unitary matrix.
% Also, S21_lead doesn't incorporate any evanescent modes.
init_data.WbyDin=3;  % W/D_{in} value > 1
[S21_lead] = generate_S21_lead(init_data);
num_modes_lead_in=size(S21_lead,2);
transmission_lead=trace(S21_lead'*S21_lead)/num_modes_lead_in;
sprintf('Transmission of S^{lead}_{21} =%f',transmission_lead)
S21_2d= ...
 S21(:,1:init_data.num_modes_prop)*S21_lead; % effective 2D transmission matrix
%------ Eigenchannel decomposition for the 2D problem ---------------------
S21_2d_prop=S21_2d(1:init_data.num_modes_prop,:);
%----- Only propagating modes are used for the eigenchannel decomposition
[U_2d,Sigma_2d,V_2d] = svd(S21_2d_prop);% Eigenchannel modes estimated 
                                        %    only for propagating modes. 
tau_2d=diag(Sigma_2d).^2;         % Transmission of eigenchannels
tau_avg_2d=mean(tau_2d);           % Average transmission
%--------------------------------------------------------------------------
c_inc=V_2d(:,1);           % Maximally transmitting (First) singular vector
c_trans=S21_2d_prop*c_inc; % Transmitted eigenchannel wavefront
sprintf('Transmission = %f',sum(abs(c_trans(1:init_data.num_modes_prop).^2)))
c_inc_for_visualisation=S21_lead*c_inc;  
                       % c_inc_for_S21 will be used for field visualisation
norm(c_inc_for_visualisation)^2
% S21_prop contains only propagating modes, hence a truncated transmission 
% matrix containing both propagating and evanescent waves on the outgoing side
% and only propagating modes on the incident side is defined for
% visualisation purposes for applying the boundary integral equation.
S21_for_eigen=S21(:,1:init_data.num_modes_prop);  
     % Incident modes only propagating, outgoing includes evanescent modes
S11_for_eigen=S11(:,1:init_data.num_modes_prop);  
      % Incident modes only propagating, outgoing includes evanescent modes

% The plot options are 'real_part','imaginary_part','magnitude'
%plot_type='real_part';  
%plot_type='imaginary_part';  
plot_type='magnitude';  
field_visualisation_and_comparison_LR(c_inc_for_visualisation,Gij_LR, ...
    S21_for_eigen,S11_for_eigen,init_data,plot_type)
side_description={'$First (Maximally $','$transmitting)$','$eigenchannel~of~S_{21}$'};
annotation('textbox', [0.005, 0.95, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on');   
%---------------- An intermediate eigenchannel ----------------------------
eig_channel_mode_no=floor(2*(num_modes_lead_in)/3);
c_inc=V_2d(:,eig_channel_mode_no);  
c_trans=S21_2d_prop*c_inc; % Transmitted eigenchannel wavefront
%c_refl=S21_2d_prop*c_inc; % Transmitted eigenchannel wavefront
sprintf('Transmission = %f',sum(abs(c_trans(1:init_data.num_modes_prop).^2)))
c_inc_for_visualisation=S21_lead*c_inc;  % c_inc_for_S21 will be used for 
                                         %             field visualisation.
norm(c_inc_for_visualisation)^2
field_visualisation_and_comparison_LR(c_inc_for_visualisation, ...
    Gij_LR,S21_for_eigen,S11_for_eigen,init_data,plot_type)
side_description={sprintf('$~%d ^{th}~eigenchannel$',eig_channel_mode_no),'$of~S_{21}$',};
annotation('textbox', [0.005, 0.95, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on'); 

%---------- Minimally transmitting eigenchannel ---------------------------
c_inc=V_2d(:,end);         % Maximally transmitting (First) singular vector
c_trans=S21_2d_prop*c_inc; % Transmitted eigenchannel wavefront
c_refl=S21_2d_prop*c_inc; % Transmitted eigenchannel wavefront
sprintf('Transmission = %f',sum(abs(c_trans(1:init_data.num_modes_prop).^2)))
c_inc_for_visualisation=S21_lead*c_inc;  % c_inc_for_S21 will be used for 
                                         %             field visualisation
norm(c_inc_for_visualisation)^2
field_visualisation_and_comparison_LR(c_inc_for_visualisation,Gij_LR, ...
    S21_for_eigen,S11_for_eigen,init_data,plot_type)
side_description={sprintf('$%d^{th} (Minimally $',size(V_2d,2)), ...
    '$transmitting)$','$eigenchannel~of~S_{21}$'};
annotation('textbox', [0.005, 0.95, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on'); 

%--------------------------------------------------------------------------
figure('position',[100 100 1500 600])
plot(1:num_modes_lead_in,tau_2d,'-*')
xlabel('$Eigenchannel~no$','Interpreter','Latex')
ylabel('$\tau$','Interpreter','Latex')
title('$Transmission~eigenchannels~of~S_{21}$','Interpreter','Latex')
set(gca,'FontSize',18);
[max(tau_2d) mean(tau_2d) min(tau_2d)]