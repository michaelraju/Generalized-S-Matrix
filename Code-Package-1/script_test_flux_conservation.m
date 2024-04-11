%------------------- Define S matrix components ---------------------------
S11_pr_pr=S11(1:init_data.num_modes_prop,1:init_data.num_modes_prop);
                  % Purely propagating S11 without evanescent modes 
S12_pr_pr=S12(1:init_data.num_modes_prop,1:init_data.num_modes_prop); 
                  % Purely propagating S21 without evanescent modes                  
S21_pr_pr=S21(1:init_data.num_modes_prop,1:init_data.num_modes_prop); 
                  % Purely propagating S21 without evanescent modes 
S22_pr_pr=S22(1:init_data.num_modes_prop,1:init_data.num_modes_prop);
                  % Purely propagating S11 without evanescent modes 
                  
S11_pr_ev=S11(1:init_data.num_modes_prop,(1+init_data.num_modes_prop):end); 
                % Incident modes evanescent, reflected modes propagating                
S12_pr_ev=S12(1:init_data.num_modes_prop,(1+init_data.num_modes_prop):end); 
                % Incident modes evanescent, transmitted modes propagating
S21_pr_ev=S21(1:init_data.num_modes_prop,(1+init_data.num_modes_prop):end); 
                % Incident modes evanescent, transmitted modes propagating
S22_pr_ev=S22(1:init_data.num_modes_prop,(1+init_data.num_modes_prop):end); 
                % Incident modes evanescent, reflected modes propagating                
                
S11_ev_ev=S11((1+init_data.num_modes_prop):end,(1+init_data.num_modes_prop):end); 
S12_ev_ev=S12((1+init_data.num_modes_prop):end,(1+init_data.num_modes_prop):end); 
S21_ev_ev=S21((1+init_data.num_modes_prop):end,(1+init_data.num_modes_prop):end); 
S22_ev_ev=S22((1+init_data.num_modes_prop):end,(1+init_data.num_modes_prop):end); 

S11_ev_pr=S11((1+init_data.num_modes_prop):end,1:init_data.num_modes_prop); 
S12_ev_pr=S12((1+init_data.num_modes_prop):end,1:init_data.num_modes_prop); 
S21_ev_pr=S21((1+init_data.num_modes_prop):end,1:init_data.num_modes_prop); 
S22_ev_pr=S22((1+init_data.num_modes_prop):end,1:init_data.num_modes_prop); 

%------------- Flux conservation : Purely propagating modes ---------------
T_pr_pr_LR=trace(S21_pr_pr*S21_pr_pr')./init_data.num_modes_prop;
R_pr_pr_LR=trace(S11_pr_pr*S11_pr_pr')./init_data.num_modes_prop;
sprintf('T_pr_pr : wave incidence left to right is %f',T_pr_pr_LR)
sprintf('R_pr_pr : wave incidence left to right is %f',R_pr_pr_LR)
sprintf('Propagating flux conservation, T_pr_pr+R_pr_pr=%f', ...
    T_pr_pr_LR+R_pr_pr_LR)


T_pr_pr_RL=trace(S12_pr_pr*S12_pr_pr')./init_data.num_modes_prop;
R_pr_pr_RL=trace(S22_pr_pr*S22_pr_pr')./init_data.num_modes_prop;
sprintf('T_pr_pr : wave incidence right to left is %f',T_pr_pr_RL)
sprintf('R_pr_pr : wave incidence right to left is %f',R_pr_pr_RL)
sprintf('Propagating flux conservation, T_pr_pr+R_pr_pr=%f', ...
    T_pr_pr_RL+R_pr_pr_RL)

figure('Position', [50 50 1200 700],'color','W');
subplot(2,2,1)
imagesc(abs(S21_pr_pr'*S21_pr_pr+S11_pr_pr'*S11_pr_pr))
title('${S^{pr,pr}_{21}}^\dagger S^{pr,pr}_{21} +{S^{pr,pr}_{11}}^\dagger S^{pr,pr}_{11}=I$',...
    'Interpreter','Latex')
colorbar
axis equal tight
set(gca,'FontSize',16)

subplot(2,2,2)
imagesc(abs(S11_pr_pr'*S12_pr_pr+S21_pr_pr'*S22_pr_pr))
title('${S^{pr,pr}_{11}}^\dagger S^{pr,pr}_{12} +{S^{pr,pr}_{21}}^\dagger S^{pr,pr}_{22}=0$',...
    'Interpreter','Latex')
colorbar
axis equal tight
set(gca,'FontSize',16)


subplot(2,2,3)
imagesc(abs(S12_pr_pr*S22_pr_pr'+S11_pr_pr*S21_pr_pr'))
title('${S^{pr,pr}_{12}}{S^{pr,pr}_{22}}^\dagger  +{S^{pr,pr}_{11}}{S^{pr,pr}_{21}}^\dagger =0$',...
    'Interpreter','Latex')
colorbar
axis equal tight
set(gca,'FontSize',16)


subplot(2,2,4)
imagesc(abs(S12_pr_pr'*S12_pr_pr+S22_pr_pr'*S22_pr_pr))
title('${S^{pr,pr}_{12}}^\dagger S^{pr,pr}_{12} +{S^{pr,pr}_{22}}^\dagger S^{pr,pr}_{22}=I$',...
    'Interpreter','Latex')
colorbar
axis equal tight
set(gca,'FontSize',16)

%------------- Evanescent portion of the S matrix  ------------------------
T_pr_ev_LR=trace(S21_pr_ev'*S21_pr_ev)./init_data.num_modes_evanes;
R_pr_ev_LR=trace(S11_pr_ev'*S11_pr_ev)./init_data.num_modes_evanes;
sprintf('T_pr_ev : wave incidence left to right is %f',T_pr_ev_LR)
sprintf('R_pr_ev : wave incidence left to right is %f',R_pr_ev_LR)
sprintf('Outgoing flux, T_pr_ev+R_pr_ev=%f',T_pr_ev_LR+R_pr_ev_LR)
opt_refl_LR=2*imag(trace(S11_ev_ev))./init_data.num_modes_evanes;
sprintf('Optical theorem part of the reflectance = %f',opt_refl_LR)
sprintf('Conservation due to evanscent wave incidence %f,%f,%f', ...
    (T_pr_ev_LR+R_pr_ev_LR),opt_refl_LR,(T_pr_ev_LR+R_pr_ev_LR)-opt_refl_LR)

T_pr_ev_RL=trace(S12_pr_ev'*S12_pr_ev)./init_data.num_modes_evanes;
R_pr_ev_RL=trace(S22_pr_ev'*S22_pr_ev)./init_data.num_modes_evanes;
sprintf('T_pr_ev : wave incidence right to left is %f',T_pr_ev_RL)
sprintf('R_pr_ev : wave incidence right to left is %f',R_pr_ev_RL)
sprintf('Outgoing flux, T_pr_ev+R_pr_ev=%f',T_pr_ev_RL+R_pr_ev_RL)
opt_refl_RL=2*imag(trace(S22_ev_ev))./init_data.num_modes_evanes;
sprintf('Optical theorem part of the reflectance = %f',opt_refl_RL)
sprintf('Conservation due to evanscent wave incidence %f,%f,%f', ...
    (T_pr_ev_RL+R_pr_ev_RL),opt_refl_RL,(T_pr_ev_RL+R_pr_ev_RL)-opt_refl_RL)


%--------------------------------------------------------------------------

figure('Position', [50 50 500 700],'color','W');
subplot(2,2,1)
imagesc(abs(S11_pr_pr'*S11_pr_ev+S21_pr_pr'*S21_pr_ev - 1i.*S11_ev_pr'))
title({'${S^{pr,pr}_{11}}^\dagger S^{pr,ev}_{11} +$', '${S^{pr,pr}_{21}}^\dagger S^{pr,ev}_{21}$','$-i {S^{ev,pr}_{11}}^\dagger $'},...
    'Interpreter','Latex')
colorbar
axis equal tight
set(gca,'FontSize',16)

subplot(2,2,2)
imagesc(abs(S11_pr_pr'*S12_pr_ev+S21_pr_pr'*S22_pr_ev - 1i.*S21_ev_pr'))
title({'${S^{pr,pr}_{11}}^\dagger S^{pr,ev}_{12} +$', '${S^{pr,pr}_{21}}^\dagger S^{pr,ev}_{22}$','$-i {S^{ev,pr}_{21}}^\dagger $'},...
    'Interpreter','Latex')
colorbar
axis equal tight
set(gca,'FontSize',16)

subplot(2,2,3)
imagesc(abs(S12_pr_pr'*S11_pr_ev+S22_pr_pr'*S21_pr_ev - 1i.*S12_ev_pr'))
title({'${S^{pr,pr}_{12}}^\dagger S^{pr,ev}_{11} +$', '${S^{pr,pr}_{22}}^\dagger S^{pr,ev}_{21}$','$-i {S^{ev,pr}_{12}}^\dagger $'},...
    'Interpreter','Latex')
colorbar
axis equal tight
set(gca,'FontSize',16)

subplot(2,2,4)
imagesc(abs(S12_pr_pr'*S12_pr_ev+S22_pr_pr'*S22_pr_ev - 1i.*S22_ev_pr'))
title({'${S^{pr,pr}_{12}}^\dagger S^{pr,ev}_{12} +$', '${S^{pr,pr}_{22}}^\dagger S^{pr,ev}_{22}$','$-i {S^{ev,pr}_{22}}^\dagger $'},...
    'Interpreter','Latex')
colorbar
axis equal tight
set(gca,'FontSize',16)
%--------------------------------------------------------------------------



figure('Position', [50 50 1200 700],'color','W');
subplot(2,2,1)
imagesc(abs(S11_pr_ev'*S11_pr_ev+...
    S21_pr_ev'*S21_pr_ev-1i.*(S11_ev_ev'-S11_ev_ev)));
title({'${S^{pr,ev}_{11}}^\dagger S^{pr,ev}_{11} + $','${S^{pr,ev}_{21}}^\dagger S^{pr,ev}_{21}$','$-i({S^{ev,ev}_{11}}^\dagger-S^{ev,ev}_{11})=0$'},...
    'Interpreter','Latex')
colorbar
axis equal tight
set(gca,'FontSize',16)

subplot(2,2,2)
imagesc(abs(S11_pr_ev'*S12_pr_ev+...
    S21_pr_ev'*S22_pr_ev-1i.*(S21_ev_ev'-S12_ev_ev)));
title({'${S^{pr,ev}_{11}}^\dagger S^{pr,ev}_{12} + $','${S^{pr,ev}_{21}}^\dagger S^{pr,ev}_{22}$','$-i({S^{ev,ev}_{21}}^\dagger-S^{ev,ev}_{12})=0$'},...
    'Interpreter','Latex')
colorbar
axis equal tight
set(gca,'FontSize',16)

subplot(2,2,3)
imagesc(abs(S12_pr_ev'*S11_pr_ev+...
    S22_pr_ev'*S21_pr_ev-1i.*(S12_ev_ev'-S21_ev_ev)));
title({'${S^{pr,ev}_{12}}^\dagger S^{pr,ev}_{11} + $','${S^{pr,ev}_{22}}^\dagger S^{pr,ev}_{21}$','$-i({S^{ev,ev}_{12}}^\dagger-S^{ev,ev}_{21})=0$'},...
    'Interpreter','Latex')
colorbar
axis equal tight
set(gca,'FontSize',16)

subplot(2,2,4)
imagesc(abs(S12_pr_ev'*S12_pr_ev+...
    S22_pr_ev'*S22_pr_ev-1i.*(S22_ev_ev'-S22_ev_ev)));
title({'${S^{pr,ev}_{12}}^\dagger S^{pr,ev}_{12} + $','${S^{pr,ev}_{22}}^\dagger S^{pr,ev}_{22}$','$-i({S^{ev,ev}_{22}}^\dagger-S^{ev,ev}_{22})=0$'},...
    'Interpreter','Latex')
colorbar
axis equal tight
set(gca,'FontSize',16)

