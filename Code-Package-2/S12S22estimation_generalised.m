function [S12,S22] = S12S22estimation_generalised(Gij_R_transpose,Gij_L_transpose,init_data)
Chi=init_data.Chi;
kz_flux=init_data.kz_flux;
kzdz=init_data.kzdz;
W=init_data.W;
jth=init_data.jth;
kydy=init_data.kydy;
dy=init_data.dy;
num_modes_prop=init_data.num_modes_prop;
num_modes=init_data.num_modes;
Nz0_R=init_data.Nz0_R;
% Nmode : scattered modes
% Mmode : incident modes
Mmode_prop=num_modes_prop;
Nmode_prop=num_modes_prop;
Mmode=num_modes;
Nmode=num_modes;
S12=((zeros(Nmode,Mmode)));
S22=((zeros(Nmode,Mmode)));

% Case 1 : Incident modes : propagating and scattered modes : propagating
Chi_n=Chi(1:Nmode_prop,W,jth.',kydy);

for mcount=1:Mmode_prop
    sprintf('Evaluating Fisher-Lee integral for m=%d',mcount)
    Chi_m=Chi(mcount,W,jth,kydy);

    prefactor_S12=2*1i*sqrt(kz_flux(mcount)).*sqrt(kz_flux(1:Nmode_prop)).*exp(-1i*kzdz(mcount)*(Nz0_R-1));
    prefactor_S22=2*1i*sqrt(kz_flux(mcount)).*sqrt(kz_flux(1:Nmode_prop)).*exp(-1i*kzdz(mcount)*(Nz0_R-1)).*exp(-1i*kzdz(1:Nmode_prop)*(Nz0_R-1));
   
    S12_temp=(Gij_L_transpose*Chi_n);
    S22_temp=(Gij_R_transpose*Chi_n);

      S12(1:Nmode_prop,mcount)=prefactor_S12.*(Chi_m*S12_temp).*(dy^2);
      S22(1:Nmode_prop,mcount)=prefactor_S22.*(Chi_m*S22_temp).*(dy^2);
      S22(1:Nmode_prop,mcount)=-(((1:Nmode_prop)==mcount).*exp(-1i*kzdz(mcount)*(Nz0_R-1)).*exp(-1i*kzdz(1:Nmode_prop)*(Nz0_R-1))).'+ S22(1:Nmode_prop,mcount);
end

% Case 2 : Incident modes : evanescent and scattered modes : evanescent 

Chi_n=Chi(1+Nmode_prop:Nmode,W,jth.',kydy);

for mcount=1+Mmode_prop:Mmode
    sprintf('Evaluating Fisher-Lee integral for m=%d',mcount)
    Chi_m=Chi(mcount,W,jth,kydy);
    
 prefactor_S12=-2*sqrt(kz_flux(mcount)).*sqrt(kz_flux(1+Nmode_prop:num_modes));
 prefactor_S22=-2*sqrt(kz_flux(mcount)).*sqrt(kz_flux(1+Nmode_prop:num_modes));
 
   
    S12_temp=(Gij_L_transpose*Chi_n);
    S22_temp=(Gij_R_transpose*Chi_n);

      S12(1+Nmode_prop:Nmode,mcount)=prefactor_S12.*(Chi_m*S12_temp).*(dy^2);
      S22(1+Nmode_prop:Nmode,mcount)=prefactor_S22.*(Chi_m*S22_temp).*(dy^2);
      S22(1+Nmode_prop:Nmode,mcount)=-((1+Nmode_prop:Nmode)==mcount).'+ S22(1+Nmode_prop:Nmode,mcount);

end

% Case 3 : Incident modes : propagating and scattered modes : evanescent
Chi_n=Chi(1+Nmode_prop:num_modes,W,jth.',kydy);

for mcount=1:Mmode_prop
    sprintf('Evaluating Fisher-Lee integral for m=%d',mcount)
    Chi_m=Chi(mcount,W,jth,kydy);
    
    prefactor_S12=2*1i*sqrt(kz_flux(mcount)).*sqrt(kz_flux(1+Nmode_prop:num_modes)).*exp(-1i*kzdz(mcount)*(Nz0_R-1));
    prefactor_S22=2*1i*sqrt(kz_flux(mcount)).*sqrt(kz_flux(1+Nmode_prop:num_modes)).*exp(-1i*kzdz(mcount)*(Nz0_R-1));
   
    S12_temp=(Gij_L_transpose*Chi_n);
    S22_temp=(Gij_R_transpose*Chi_n);

      S12(1+Nmode_prop:num_modes,mcount)=prefactor_S12.*(Chi_m*S12_temp).*(dy^2);
      S22(1+Nmode_prop:num_modes,mcount)=prefactor_S22.*(Chi_m*S22_temp).*(dy^2);
end

%Case 4 : Incident modes : evanescent and scattered modes : propagating 
Chi_n=Chi(1:Nmode_prop,W,jth.',kydy);

for mcount=1+Mmode_prop:Mmode
    sprintf('Evaluating Fisher-Lee integral for m=%d',mcount)
    Chi_m=Chi(mcount,W,jth,kydy);
    
 prefactor_S12=-2*sqrt(kz_flux(mcount)).*sqrt(kz_flux(1:Nmode_prop));
 prefactor_S22=-2*sqrt(kz_flux(mcount)).*sqrt(kz_flux(1:Nmode_prop)).*exp(-1i*kzdz(1:Nmode_prop)*(Nz0_R-1));
 
   
    S12_temp=(Gij_L_transpose*Chi_n);
    S22_temp=(Gij_R_transpose*Chi_n);

      S12(1:Nmode_prop,mcount)=prefactor_S12.*(Chi_m*S12_temp).*(dy^2);
      S22(1:Nmode_prop,mcount)=prefactor_S22.*(Chi_m*S22_temp).*(dy^2);
end

%------------------------- Plotting ---------------------------------------
num_modes_prop=init_data.num_modes_prop;
num_modes=init_data.num_modes;
num_modes_evanes=init_data.num_modes_evanes;
ticks= [1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop num_modes_evanes];

figure('Position', [0 0 1600 800],'color','W');

subplot(2,2,1)
colormap jet
imagesc(abs(S12).^2)
axis equal tight
colorbar
xlabel('$incident~mode~no$','Interpreter','Latex')
ylabel('$transmitted~mode~no$','Interpreter','Latex')
xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop num_modes]);
yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop num_modes]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
caxis([min(min((abs(S12).^2))) max(max((abs(S12).^2)))])
set(gca,'ColorScale','log')
title('$|S_{12}|^2$','Interpreter','Latex')
set(gca,'FontSize',18) 

subplot(2,2,2)
colormap jet
%imagesc(angle(S12).*((abs(S12).^2)> 10^-10)) % set a threshold of 10^(-10) for 
                                      % on the magnitude squared to avoid plotting  
                                      % phase for magnitude values less
                                      % than the threshold.
imagesc(angle(S12))                                      
axis equal tight
colorbar
caxis([-pi pi])
xlabel('$incident~mode~no$','Interpreter','Latex')
ylabel('$transmitted~mode~no$','Interpreter','Latex')
xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop num_modes]);
yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop num_modes]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
title('$Phase(S_{12})$','Interpreter','Latex')
set(gca,'FontSize',18) 


subplot(2,2,3)
colormap jet
imagesc(abs(S22).^2)
axis equal tight
colorbar
xlabel('$incident~mode~no$','Interpreter','Latex')
ylabel('$reflected~mode~no$','Interpreter','Latex')
xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop num_modes]);
yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop num_modes]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
caxis([min(min((abs(S22).^2))) max(max((abs(S22).^2)))])
set(gca,'ColorScale','log')
title('$|S_{22}|^2$','Interpreter','Latex')
set(gca,'FontSize',18) 

subplot(2,2,4)
colormap jet
%imagesc(angle(S22).*((abs(S22).^2)>10^-10))% set a threshold of 10^(-10) for 
                                      % on the magnitude to avoid plotting  
                                      % phase for magnitude values less
                                      % than the threshold.
imagesc(angle(S22))                                      
axis equal tight
caxis([-pi pi])
colorbar
xlabel('$incident~mode~no$','Interpreter','Latex')
ylabel('$reflected~mode~no$','Interpreter','Latex')
xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop num_modes]);
yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop num_modes]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
title('$Phase(S_{22})$','Interpreter','Latex')
set(gca,'FontSize',18) 
%---------------------- Annotating the figure ------------------------------
side_description={sprintf('No of modes'),sprintf('propagating = %d',init_data.num_modes_prop), ...
                  sprintf('evanescent= %d',init_data.num_modes_evanes)};
annotation('textbox', [0.005, 0.9, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on');   
end

