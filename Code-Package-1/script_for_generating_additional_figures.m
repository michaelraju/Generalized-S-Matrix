%------------------- Generate Fig 1 for the paper -------------------------
FontSizeVal=22;
%------ plotting a propagating eigenmode and two evanescent modes ---------
mcount_prop=3;
prop_mode=init_data.find_prop_modes_LR(mcount_prop,init_data.kzdz, ...
    init_data.dz,init_data.Nmat,init_data.W,init_data.Mmat,init_data.kydy);

mcount_evan=init_data.num_modes_prop+2; % evanescent mode number

evan_mode_L=init_data.find_evanes_modes(mcount_evan,init_data.kzdz, ...
    init_data.dz,init_data.Nmat,init_data.W,init_data.Mmat, ...
    init_data.kydy,init_data.Nz0_L);
evan_mode_R=init_data.find_evanes_modes(mcount_evan,init_data.kzdz,...
    init_data.dz,init_data.Nmat,init_data.W,init_data.Mmat,init_data.kydy,...
    init_data.Nz0_R);


figure('position',[100 100 1800 800])
subplot(1,5,1)
imagesc([0 init_data.kref*init_data.dz*(init_data.Nz-1)],[0 ...
    init_data.kref*init_data.dy*(init_data.Ny-1)],real(prop_mode))
axis xy equal tight
colormap wavecolormap
colorbar 
xlabel('$k_{ref}z$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
title('$real[\psi^{+}_{m,pr}(z,y)]$','Interpreter','Latex')
set(gca,'FontSize',FontSizeVal) 

subplot(1,5,2)
imagesc([0 init_data.kref*init_data.dz*(init_data.Nz-1)],[0 ...
    init_data.kref*init_data.dy*(init_data.Ny-1)],real(evan_mode_L))
axis xy equal tight
colormap wavecolormap
colorbar 
xlabel('$k_{ref}z$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
title('$\psi^{+}_{m,ev}(z,y,z_L)$','Interpreter','Latex')
set(gca,'FontSize',FontSizeVal) 

subplot(1,5,3)
imagesc([0 init_data.kref*init_data.dz*(init_data.Nz-1)],[0 ...
    init_data.kref*init_data.dy*(init_data.Ny-1)],real(evan_mode_R))
axis xy equal tight
colormap wavecolormap
colorbar 
xlabel('$k_{ref}z$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
title('$\psi^{-}_{m,ev}(z,y,z_R)$','Interpreter','Latex')
set(gca,'FontSize',FontSizeVal) 

subplot(1,5,4)
Ny_loc=floor(init_data.Ny/2)+2;  % Choose a location 
[(Ny_loc)*init_data.kref*init_data.dy]
imagesc([0 init_data.kref*init_data.dz*(init_data.Nz-1)],[0 ...
    init_data.kref*init_data.dy*(init_data.Ny-1)],real(G0ij_LR(:,:,Ny_loc)))
axis xy equal tight
caxis([-max(max(real(G0ij_LR(:,:,Ny_loc)))) max(max(real(G0ij_LR(:,:,Ny_loc))))])
colormap wavecolormap
colorbar 
xlabel('$k_{ref}z$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
title('$real[G_0(r,r_0)]$','Interpreter','Latex')
set(gca,'FontSize',FontSizeVal) 

subplot(1,5,5)
G_for_plotting=reshape(Gij_LR(:,:,Ny_loc), ...
    init_data.Ny,init_data.Nz);
imagesc([0 init_data.kref*init_data.dz*(init_data.Nz-1)],[0 ...
    init_data.kref*init_data.dy*(init_data.Ny-1)],real(G_for_plotting))
axis xy equal tight
caxis([-max(real(G_for_plotting(:))) max(real(G_for_plotting(:)))])
colormap wavecolormap
colorbar 
xlabel('$k_{ref}z$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
title('$real[G(r,r_0)]$','Interpreter','Latex')
set(gca,'FontSize',FontSizeVal) 
%--------------------- Annotation -----------------------------------------
side_description={'Propagating',...
                  sprintf('$mode~no~%d$',mcount_prop), ...
                  'Evanescent',...
                  sprintf('$mode~no~%d$',mcount_evan-init_data.num_modes_prop),...
                  'No~of~evanescent',...
                  sprintf('modes~taken=~%d',init_data.num_modes_evanes)};
annotation('textbox', [0.005, 0.95, 0.001, 0.001], 'string', ...
   side_description,'FontSize',22, ...
   'Interpreter','latex','FitBoxToText','on');   


%------------------ Plot the disorder -------------------------------------
estimate_disorder_spatial_correlation(eps_profile,init_data);


%---- save the figure in vector format ------------------------------------
%exportgraphics(gcf,'eigenmodes_and_greens_fun.eps','ContentType','vector')



%------------visualising modes for finite beam incidence ------------------
% m_count=3;
% c_inc=zeros(size(S21_lead,2),1);
% c_inc(m_count)=1;
% c_inc=c_inc./norm(c_inc);
% c_inc_for_visualisation=S21_lead*c_inc;     
% sprintf('Lead Transmission = %f',sum(abs(c_inc_for_visualisation(1:init_data.num_modes_prop).^2)))
% plot_type='magnitude';  
% field_visualisation_and_comparison_LR(c_inc_for_visualisation,Gij_LR,S21_for_eigen,S11_for_eigen,init_data,plot_type)
% title('$|Total~field|$','Interpreter','Latex');
% set(gca,'FontSize',FontSizeVal);
% plot_type='real_part';  
% %--- a quick plotting for free space propagation without any disorder -----
% % Taking transmission matrix as an identity matrix and reflection matrix to
% % be zero for quick visualisation of waves in a disorder-less region
% field_visualisation_and_comparison_LR(c_inc_for_visualisation,G0ij_LR,eye(size(S21_for_eigen)),zeros(size(S11_for_eigen)),init_data,plot_type)
% set(gca,'FontSize',FontSizeVal);