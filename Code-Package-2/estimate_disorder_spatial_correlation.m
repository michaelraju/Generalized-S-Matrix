function [kreflc_main_numerical,var_disorder] = estimate_disorder_spatial_correlation(eps_profile,init_data)
FontSizeVal=18;
no_of_bins=50;
Msample=init_data.Ny-2;
Nsample=init_data.Nz-2*init_data.offset;
eps_profile_non_zero=reshape(eps_profile(find(eps_profile)),Msample,Nsample);
sample_corr=abs(fftshift(ifft2(fft2(eps_profile_non_zero).*conj(fft2(eps_profile_non_zero)))))./(Msample*Nsample);
[Mcentre,Ncentre]=ind2sub([Msample Nsample],find(sample_corr==max(sample_corr(:))));
samp_main_corr=sample_corr(Mcentre,Ncentre:end);
%-------------- Interpolate for the 1/e correlation value -----------------
one_by_e_value=samp_main_corr(1).*1/exp(1);
ind_before=find(samp_main_corr>one_by_e_value,1,'last'); 
ind_after=find(samp_main_corr<one_by_e_value,1,'first');
x_indices=0:length(samp_main_corr)-1;

N_lc_main_numerical = ...
    interp1([samp_main_corr(ind_before) samp_main_corr(ind_after)], ...
    [x_indices(ind_before) x_indices(ind_after)], ...
    samp_main_corr(1).*1/exp(1));

%--------------Estimate the dimensionless spatial correlation length ------ 
kreflc_main_numerical=init_data.kref*init_data.dz*N_lc_main_numerical;
%----------------- Plot the disorder properties ---------------------------
figure('Position', [0 0 1100 700],'color','W');
subplot(2,2,4)
plot((0:length(samp_main_corr)-1).*(init_data.kref*init_data.dz), ...
    samp_main_corr,'-*','MarkerSize',8,'LineWidth',1)
xlabel('$k_{ref}|z-z_0|$','Interpreter','Latex')
ylabel('$R(k_{ref}|z-z_0|)$','Interpreter','Latex')
set(gca,'FontSize',FontSizeVal); 

subplot(2,2,[1 3])
colormap(jet);
imagesc([0 init_data.kref*init_data.dz*(init_data.Nz-1)], ...
    [0 init_data.kref*init_data.dy*(init_data.Ny-1)], ...
    (eps_profile)); 
xlabel('$k_{ref}z$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
axis xy 
axis equal tight
set(gca,'FontSize',FontSizeVal) 
colorbar
sprintf('No of particles is %d',length(find(eps_profile)))
title('$\delta \epsilon_r(z,y)$','Interpreter','Latex')
non_zero_eps=eps_profile(find(eps_profile));
[hist_data_main,x_hist] = histcounts(non_zero_eps,no_of_bins, ...
    'Normalization', 'pdf');
mid_point = x_hist(1:end-1)+mean(diff(x_hist))/2;   
dtau=(mid_point(2)-mid_point(1));
[sum(hist_data_main.*dtau) 1]  % Histogram
[var(non_zero_eps)  ((mid_point(end)-mid_point(1))^2)/12 samp_main_corr(1)]

subplot(2,2,2)
plot(mid_point,hist_data_main,'-*','MarkerSize',8,'LineWidth',1);
xlabel('$\delta\epsilon_r$','Interpreter','Latex')
ylabel('$p(\delta\epsilon_r)$','Interpreter','Latex')
title('$P.D.F~of~\delta \epsilon(r)$','Interpreter','Latex')
set(gca,'FontSize',FontSizeVal) 
%-------------------- Annotation ------------------------------------------
var_disorder=var(non_zero_eps); 
side_description={'$var(disorder)$', ...
                  sprintf('$=%.4f$',var_disorder), ...
                  '$k_{ref}l_c=$',...
                  sprintf('$=%0.4f$',kreflc_main_numerical), ...
                  '$R_e={(k_{ref}l_c)}^2 var_{disorder}$',...
                  sprintf('=%.3f',var_disorder*kreflc_main_numerical^2)};
annotation('textbox', [0.005, 0.95, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on');   
%--------------------------------------------------------------------------
end

