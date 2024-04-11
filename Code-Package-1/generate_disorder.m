function eps_profile = generate_disorder(init_data)
Ny=init_data.Ny;
Nz=init_data.Nz;
sigma_disorder=init_data.sigma_disorder;
n0=init_data.n0;
disorder_corr_flag=init_data.disorder_corr_flag;
%---------------- Obtain the dielectric constant spatial profile ---------- 
eps_profile=zeros(Ny,Nz);
rng('shuffle');          %  random number seed based on time 
offset=init_data.offset; % offset from sides of the boundary where 
                         % particles are present

if disorder_corr_flag==1   % For spatially correlated disorder
 if init_data.disorder_perturbation_strength_parameter==0
 sprintf('For correlated disorders, disorder_perturbation_strength_parameter has to be non-zero !!! ')
 pause ;
 end
[corr_disorder_main_uniform]=correlated_disorder_fun(init_data);
eps_profile(:,offset+1:end-offset)=corr_disorder_main_uniform(:,offset+1:end-offset);
else                       % For spatially uncorrelated disorder
eps_profile(:,offset+1:end-offset)=sigma_disorder.*(2.*rand(Ny,Nz-2*offset)-1); 
end

eps_profile(1,offset+1:end-1)=0; % Force the transverse boundary dielectric
                                 % perturbation to be zero (because its
                                 % the grid boundary)
eps_profile(end,offset+1:end-1)=0; % Force the transverse boundary dielectric 
                                   % perturbation to be zero (because its 
                                   % the grid boundary) 
if min(min(eps_profile + n0^2)) < 1
    disp('Warning !! relative dielectric constant less than 1');
    pause;
end

%save('eps_profile.mat','eps_profile');
%---------------------------- plot ----------------------------------------
% figure('Position', [0 0 1000 1000],'color','W');
% subplot(1,2,1)
% colormap(jet);
% imagesc([0 kref*dz*(Nz-1)],[0 kref*dy*(Ny-1)],(eps_profile));  % Real part of eigenmode
% xlabel('$k_{ref}z$','Interpreter','Latex')
% ylabel('$k_{ref}y$','Interpreter','Latex')
% axis xy 
% axis equal tight
% set(gca,'FontSize',18) 
% colorbar
% sprintf('No of particles is %d',length(find(eps_profile)))
% title('$\Delta \epsilon(r)$','Interpreter','Latex')
% %-------- Estimate disorder probability distribution ----------------------
% no_of_bins=50;
% non_zero_eps=eps_profile(find(eps_profile));
% [hist_data_main,x_hist] = histcounts(non_zero_eps,no_of_bins,'Normalization', 'pdf');
% mid_point = x_hist(1:end-1)+mean(diff(x_hist))/2;   
% dtau=(mid_point(2)-mid_point(1));
% [sum(hist_data_main.*dtau) 1]  % Histogram
% %--------------------------------------------------------------------------
% subplot(1,2,2)
% plot(mid_point,hist_data_main,'-*','MarkerSize',3);
% xlabel('$\delta\epsilon_r$','Interpreter','Latex')
% ylabel('$p(\delta\epsilon_r)$','Interpreter','Latex')
% title('$P.D.F~of~\Delta \epsilon(r)$','Interpreter','Latex')
% set(gca,'FontSize',18) 
%--------------------------------------------------------------------------
end

