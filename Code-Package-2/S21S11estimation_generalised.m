function [S21,S11] = S21S11estimation_generalised(Gij_L_transpose,Gij_R_transpose,init_data)
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
%%%% Nmode : output modes
%%%% Mmode : input modes
Mmode_prop=num_modes_prop;
Nmode_prop=num_modes_prop;
Mmode=num_modes;
Nmode=num_modes;
S21=((zeros(Nmode,Mmode)));
S11=((zeros(Nmode,Mmode)));

% Case 1 : Incident modes : propagating and scattered modes : propagating
Chi_n=Chi(1:Nmode_prop,W,jth.',kydy);

for mcount=1:Mmode_prop
    sprintf('Evaluating Fisher-Lee integral for m=%d',mcount)
    Chi_m=Chi(mcount,W,jth,kydy);

    prefactor_S21=2*1i*sqrt(kz_flux(mcount)).*sqrt(kz_flux(1:Nmode_prop)).*exp(-1i*kzdz(1:Nmode_prop)*(Nz0_R-1));
    prefactor_S11=2*1i*sqrt(kz_flux(mcount)).*sqrt(kz_flux(1:Nmode_prop));
   
    S21_temp=(Gij_R_transpose*Chi_n);
    S11_temp=(Gij_L_transpose*Chi_n);

      S21(1:Nmode_prop,mcount)=prefactor_S21.*(Chi_m*S21_temp).*(dy^2);
      S11(1:Nmode_prop,mcount)=prefactor_S11.*(Chi_m*S11_temp).*(dy^2);
      S11(1:Nmode_prop,mcount)=-((1:Nmode_prop)==mcount).'+ S11(1:Nmode_prop,mcount);
end

% Case 2 : Incident modes : evanescent and scattered modes : evanescent 

Chi_n=Chi(1+Nmode_prop:Nmode,W,jth.',kydy);

for mcount=1+Mmode_prop:Mmode
    sprintf('Evaluating Fisher-Lee integral for m=%d',mcount)
    Chi_m=Chi(mcount,W,jth,kydy);
    
 prefactor_S21=-2*sqrt(kz_flux(mcount)).*sqrt(kz_flux(1+Nmode_prop:num_modes));
 prefactor_S11=-2*sqrt(kz_flux(mcount)).*sqrt(kz_flux(1+Nmode_prop:num_modes));
 
   
    S21_temp=(Gij_R_transpose*Chi_n);
    S11_temp=(Gij_L_transpose*Chi_n);

      S21(1+Nmode_prop:Nmode,mcount)=prefactor_S21.*(Chi_m*S21_temp).*(dy^2);
      S11(1+Nmode_prop:Nmode,mcount)=prefactor_S11.*(Chi_m*S11_temp).*(dy^2);
      S11(1+Nmode_prop:Nmode,mcount)=-((1+Nmode_prop:Nmode)==mcount).'+ S11(1+Nmode_prop:Nmode,mcount);

end

% Case 3 : Incident modes : propagating and scattered modes : evanescent
Chi_n=Chi(1+Nmode_prop:num_modes,W,jth.',kydy);

for mcount=1:Mmode_prop
    sprintf('Evaluating Fisher-Lee integral for m=%d',mcount)
    Chi_m=Chi(mcount,W,jth,kydy);
    
    prefactor_S21=2*1i*sqrt(kz_flux(mcount)).*sqrt(kz_flux(1+Nmode_prop:num_modes));
    prefactor_S11=2*1i*sqrt(kz_flux(mcount)).*sqrt(kz_flux(1+Nmode_prop:num_modes));
   
    S21_temp=(Gij_R_transpose*Chi_n);
    S11_temp=(Gij_L_transpose*Chi_n);

      S21(1+Nmode_prop:num_modes,mcount)=prefactor_S21.*(Chi_m*S21_temp).*(dy^2);
      S11(1+Nmode_prop:num_modes,mcount)=prefactor_S11.*(Chi_m*S11_temp).*(dy^2);
end

%Case 4 : Incident modes : evanescent and scattered modes : propagating 
Chi_n=Chi(1:Nmode_prop,W,jth.',kydy);

for mcount=1+Mmode_prop:Mmode
    sprintf('Evaluating Fisher-Lee integral for m=%d',mcount)
    Chi_m=Chi(mcount,W,jth,kydy);
    
 prefactor_S21=-2*sqrt(kz_flux(mcount)).*sqrt(kz_flux(1:Nmode_prop)).*exp(-1i*kzdz(1:Nmode_prop)*(Nz0_R-1));
 prefactor_S11=-2*sqrt(kz_flux(mcount)).*sqrt(kz_flux(1:Nmode_prop));
 
   
    S21_temp=(Gij_R_transpose*Chi_n);
    S11_temp=(Gij_L_transpose*Chi_n);

      S21(1:Nmode_prop,mcount)=prefactor_S21.*(Chi_m*S21_temp).*(dy^2);
      S11(1:Nmode_prop,mcount)=prefactor_S11.*(Chi_m*S11_temp).*(dy^2);
end

%------------------------- Plotting ---------------------------------------
num_modes_prop=init_data.num_modes_prop;
num_modes=init_data.num_modes;
num_modes_evanes=init_data.num_modes_evanes;
ticks= [1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop num_modes_evanes];

figure('Position', [0 0 1600 800],'color','W');
subplot(2,2,1)
colormap jet
imagesc(abs(S21).^2)
axis equal tight
colorbar
xlabel('$incident~mode~no$','Interpreter','Latex')
ylabel('$transmitted~mode~no$','Interpreter','Latex')
xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop num_modes]);
yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop num_modes]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
caxis([min(min((abs(S21).^2))) max(max((abs(S21).^2)))])
set(gca,'ColorScale','log')
title('$|S_{21}|^2$','Interpreter','Latex')
set(gca,'FontSize',18) 

subplot(2,2,2)
colormap jet
%imagesc(angle(S21).*((abs(S21).^2)> 10^-10)) % set a threshold of 10^(-10) for 
                                      % on the magnitude squared to avoid plotting  
                                      % phase for magnitude values less
                                      % than the threshold.
imagesc(angle(S21))                                      
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
title('$Phase(S_{21})$','Interpreter','Latex')
set(gca,'FontSize',18) 


subplot(2,2,3)
colormap jet
imagesc((abs(S11).^2))                                     
axis equal tight
colorbar
xlabel('$incident~mode~no$','Interpreter','Latex')
ylabel('$reflected~mode~no$','Interpreter','Latex')
xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop num_modes]);
yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop num_modes]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
caxis([min(min((abs(S11).^2))) max(max((abs(S11).^2)))])
set(gca,'ColorScale','log')
title('$|S_{11}|^2$','Interpreter','Latex')
set(gca,'FontSize',18) 

subplot(2,2,4)
colormap jet
%imagesc(angle(S11).*((abs(S11).^2)>10^-10))% set a threshold of 10^(-10) for 
                                      % on the magnitude to avoid plotting  
                                      % phase for magnitude values less
                                      % than the threshold.
imagesc(angle(S11))                                      
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
set(gca,'FontSize',18) 
title('$Phase(S_{11})$','Interpreter','Latex')
%---------------------- Annotating the figure ------------------------------
side_description={sprintf('No of modes'),sprintf('propagating = %d',init_data.num_modes_prop), ...
                  sprintf('evanescent= %d',init_data.num_modes_evanes)};
annotation('textbox', [0.005, 0.9, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on');   
end

