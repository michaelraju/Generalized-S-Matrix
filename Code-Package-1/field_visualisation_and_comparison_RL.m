function field_visualisation_and_comparison_RL(c_inc,Gij_RL,S12,S22,init_data,plot_type)
Nz0_L=init_data.Nz0_L;  % z index for the left boundary
Nz0_R=init_data.Nz0_R;  % z index for the right boundary
%------------------------- Testing BI equation ----------------------------
c_trans=S12*c_inc;
c_refl=S22*c_inc;
evanes_flag=~(isempty(find(c_inc(1+init_data.num_modes_prop:end),1)));
if evanes_flag==0
T=norm(c_trans(1:init_data.num_modes_prop))^2;
R=norm(c_refl(1:init_data.num_modes_prop))^2;
sprintf('Transmission=%f, Reflection=%f, R+T=%f',T,R,R+T)
end
E_trans=zeros(size(init_data.Z(:,end))); % Transmitted field
E_refl=zeros(size(init_data.Z(:,end)));  % Reflected field
E_inc=zeros(size(init_data.Z));          % Incident field

%------------- Estimate the incident wave ---------------------------------
for mcount=1:init_data.num_modes_prop % The propagating modes
E_inc = E_inc+c_inc(mcount).*init_data.find_prop_modes_RL(mcount, ...
    init_data.kzdz,init_data.dz,init_data.Nmat,init_data.W,init_data.Mmat,...
    init_data.kydy);
end

if length(c_inc)> init_data.num_modes_prop % The evanescent modes
    for mcount=1+init_data.num_modes_prop:init_data.num_modes
    E_inc = E_inc+c_inc(mcount).*init_data.find_evanes_modes(mcount,...
        init_data.kzdz,init_data.dz,init_data.Nmat,init_data.W,init_data.Mmat,...
        init_data.kydy,init_data.Nz0_R);
    end
end
%--------------------- Estimate reflected and transmitted waves -----------
for mcount=1:init_data.num_modes_prop
E_trans = E_trans+c_trans(mcount).*init_data.find_prop_modes_RL(mcount,...
    init_data.kzdz,init_data.dz,init_data.Nmat(:,Nz0_L),init_data.W,...
    init_data.Mmat(:,Nz0_L),init_data.kydy);
E_refl = E_refl+c_refl(mcount).*init_data.find_prop_modes_LR(mcount,...
    init_data.kzdz,init_data.dz,init_data.Nmat(:,Nz0_R),init_data.W,...
    init_data.Mmat(:,Nz0_R),init_data.kydy);
end

if length(c_trans)> init_data.num_modes_prop
   for mcount=1+init_data.num_modes_prop:init_data.num_modes
   E_trans = E_trans+c_trans(mcount).*init_data.find_evanes_modes(mcount,...
       init_data.kzdz,init_data.dz,init_data.Nmat(:,Nz0_L),init_data.W,...
       init_data.Mmat(:,Nz0_L),init_data.kydy,init_data.Nz0_L);
   E_refl = E_refl+c_refl(mcount).*init_data.find_evanes_modes(mcount,...
       init_data.kzdz,init_data.dz,init_data.Nmat(:,Nz0_R),init_data.W,...
       init_data.Mmat(:,Nz0_R),init_data.kydy,init_data.Nz0_R);
   end
end
%-------- testing S12 and Boundary integral equation ----------------------
FontSizeVal=18;
figure('Position', [100 100 1620 780],'color','W');
subplot(1,8,7)
plot(abs(E_trans),init_data.kref*init_data.dy*(init_data.jth-1),'-*b');
hold on
total_field=find_total_field_RL(Gij_RL,c_inc,init_data);
plot(abs(total_field(:,1)),init_data.kref*init_data.dy*(init_data.jth-1),'*r');
ylabel('$k_{ref}y$','Interpreter','Latex')
title('$|\tilde{E}_{trans}|$','Interpreter','Latex')
legend('$Using~S_{12}$','$B.I~method$','Interpreter','Latex');
set(gca,'FontSize',FontSizeVal)

subplot(1,8,8)
plot(abs(E_refl),init_data.kref*init_data.dy*(init_data.jth-1),'-*b');
hold on
plot(abs(total_field(:,end)-E_inc(:,end)),...
    init_data.kref*init_data.dy*(init_data.jth-1),'*r');
ylabel('$k_{ref}y$','Interpreter','Latex')
title({'$|\tilde{E}_{refl}|=$','$|\tilde{E}_{total}-\tilde{E}_{inc}|$'},...
    'Interpreter','Latex')
legend('$Using~S_{22}$','$B.I~method$','Interpreter','Latex');
set(gca,'FontSize',FontSizeVal)
%--------------------------------------------------------------------------
subplot(1,8,[1 2 3])
if strcmp(plot_type,'real_part')==1
imagesc([0 init_data.kref*init_data.dz*(init_data.Nz-1)],[0 ...
    init_data.kref*init_data.dy*(init_data.Ny-1)],real(total_field)) 
                                                   % Real part of eigenmode
colormap wavecolormap
title('$real(Total~field,~\tilde{E})$','Interpreter','Latex')
elseif strcmp(plot_type,'imaginary_part')==1
imagesc([0 init_data.kref*init_data.dz*(init_data.Nz-1)],[0 ...
    init_data.kref*init_data.dy*(init_data.Ny-1)],imag(total_field)) 
                                                   % imag part of eigenmode
colormap wavecolormap
title('$imag(Total~field,~\tilde{E})$','Interpreter','Latex')
elseif strcmp(plot_type,'magnitude')==1
imagesc([0 init_data.kref*init_data.dz*(init_data.Nz-1)],[0 ...
    init_data.kref*init_data.dy*(init_data.Ny-1)],abs(total_field)) 
                                              % Magnitude part of eigenmode
colormap hot
title('$|Total~field,~\tilde{E}|$','Interpreter','Latex')
end
start_slab_nz=init_data.Nz0_L+init_data.offset;  
end_slab_nz=init_data.Nz0_R-init_data.offset;  
hold on
line([(start_slab_nz-1)*init_data.kref*init_data.dz ...
    (end_slab_nz-1)*init_data.kref*init_data.dz  ...
    (end_slab_nz-1)*init_data.kref*init_data.dz ...
    (start_slab_nz-1)*init_data.kref*init_data.dz  ...
    (start_slab_nz-1)*init_data.kref*init_data.dz ],[0 0 ...
    (init_data.Ny-1)*init_data.kref*init_data.dy ...
    (init_data.Ny-1)*init_data.kref*init_data.dy  0],...
    'LineWidth',2,'color','white')
axis xy equal tight
xlabel('$k_{ref}z$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
set(gca,'FontSize',FontSizeVal) 
colorbar
%-------------------------- Plot the TFSF ---------------------------------
TFSF_field=total_field;
TFSF_field(:,end_slab_nz:end)=TFSF_field(:,end_slab_nz:end)- ...
    E_inc(:,end_slab_nz:end);

subplot(1,8,[4 5 6])
if strcmp(plot_type,'real_part')==1
imagesc([0 init_data.kref*init_data.dz*(init_data.Nz-1)],[0 ...
    init_data.kref*init_data.dy*(init_data.Ny-1)],...
    real(TFSF_field)) % Real part of eigenmode
colormap wavecolormap
title('$real(Total~field/Scatter~Field,~\tilde{E})$','Interpreter','Latex')
elseif strcmp(plot_type,'imaginary_part')==1
imagesc([0 init_data.kref*init_data.dz*(init_data.Nz-1)],[0 ...
    init_data.kref*init_data.dy*(init_data.Ny-1)],imag(TFSF_field)) 
                                                   % imag part of eigenmode
colormap wavecolormap
title('$imag(Total~field/Scatter~Field,~\tilde{E})$','Interpreter','Latex')
elseif strcmp(plot_type,'magnitude')==1
imagesc([0 init_data.kref*init_data.dz*(init_data.Nz-1)],[0 ...
    init_data.kref*init_data.dy*(init_data.Ny-1)],abs(TFSF_field)) 
                                              % magnitude part of eigenmode
colormap hot
title('$|Total~field/Scatter~Field,~\tilde{E}|$','Interpreter','Latex')
end
start_slab_nz=init_data.Nz0_L+init_data.offset;  
end_slab_nz=init_data.Nz0_R-init_data.offset;  
hold on
line([(end_slab_nz-1)*init_data.kref*init_data.dz ...
    (start_slab_nz-1)*init_data.kref*init_data.dz ...
    (start_slab_nz-1)*init_data.kref*init_data.dz ...
    (end_slab_nz-1)*init_data.kref*init_data.dz],[...
    (init_data.Ny-1)*init_data.kref*init_data.dy  ...
    (init_data.Ny-1)*init_data.kref*init_data.dy  0 0],...
    'LineWidth',1,'color','white')
hold on
line([(end_slab_nz-1)*init_data.kref*init_data.dz ...
    (end_slab_nz-1)*init_data.kref*init_data.dz ],[...
    0 (init_data.Ny-1)*init_data.kref*init_data.dy],...
    'LineWidth',1,'color','Green','LineStyle','--')
axis xy equal tight
xlabel('$k_{ref}z$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
set(gca,'FontSize',FontSizeVal) 
colorbar

%---------------------- create plot annotations ---------------------------
if evanes_flag==0
annotation('textbox', [0.005, 0.8, 0.001, 0.001], 'string', ...
    {'$Wave~incidence~:$','$from~right~to~left$',sprintf('$T=%.4f$',T),...
    sprintf('$R=%.4f$',R),sprintf('$R+T=%.4f$',R+T) ...
    },'FontSize',FontSizeVal,'Interpreter','Latex','FitBoxToText','on');
annotation('textbox', [0.005, 0.6, 0.001, 0.001], 'string', ...
   'TF/SF boundary \color[rgb]{0.0, 1.0, 0.0} ---','FontSize',FontSizeVal-2, ...
   'Interpreter','tex','FitBoxToText','on');
else 
annotation('textbox', [0.005, 0.8, 0.001, 0.001], 'string', ...
    {'$Wave~incidence~:$','$from~right~to~left$'},'FontSize',FontSizeVal, ...
    'Interpreter','Latex','FitBoxToText','on');
annotation('textbox', [0.005, 0.6, 0.001, 0.001], 'string', ...
   'TF/SF boundary \color[rgb]{0.0, 1.0, 0.0} ---','FontSize',FontSizeVal-2, ...
   'Interpreter','tex','FitBoxToText','on');   
end

end

