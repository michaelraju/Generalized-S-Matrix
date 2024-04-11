FontSizeVal=18;
%------ Focussing of a propagating mode focussing via phase conjugation --
S12_pr_pr=S12(1:init_data.num_modes_prop,1:init_data.num_modes_prop); 
                  % Purely propagating S21 without evanescent modes                  
S21_pr_pr=S21(1:init_data.num_modes_prop,1:init_data.num_modes_prop); 
                  % Purely propagating S21 without evanescent modes 

mode_num_to_focus_pr=6;                   % The ith mode to be focussed
I_m_pr=zeros(init_data.num_modes_prop,1); % Column vector to choose the mode 
                                       %     to focus
I_m_pr(mode_num_to_focus_pr)=1;
%--------- Start by exciting the mode to be focussed from the right side
% of the disorder, to be phase conjugated again from the left to focus back
% on the right side. 
c_L_pr_pr=S12_pr_pr*I_m_pr; % outgoing wave on the left due the single mode 
                      % incidence from the right
c1=1/norm(c_L_pr_pr);  % Flux normalization term for phase conjugation to the
                    % right side. 
c_R_opt_pr_pr=S21_pr_pr*(c1.*(conj(S12_pr_pr))*I_m_pr); 
                     %Phase conjugation to right side of the slab
                     % Here, (c1.*(conj(S12_pr_pr))*I_m_pr) is the incident
                     % flux-normalized phase conjugate wave from the left.  

T_focus_pr_pr=abs(c_R_opt_pr_pr(mode_num_to_focus_pr))^2;
                     % Transmission just due to the focus
T_background_pr_pr=sum(abs(c_R_opt_pr_pr(1:end ~= mode_num_to_focus_pr)).^2);                    

figure('Position', [50 50 900 500],'color','W');
plot(1:init_data.num_modes_prop,abs(c_R_opt_pr_pr).^2,'-*')
xlabel('$Mode~number$','Interpreter','Latex')
ylabel('Transmission')
title({'Optimal~transmission through a single disorder : ',...
    'Single~mode~focussing~of~a~\textbf{propagating}~mode'},...
    'Interpreter','Latex')
annotation('textbox', [0.3, 0.8, 0.001, 0.001], 'string', ...
    {sprintf('$Prop~mode~num~focussed=%d$',mode_num_to_focus_pr),...
    sprintf('$T_{focus} = %.3f$',T_focus_pr_pr), ...
    sprintf('$T_{background} = %.3f$',T_background_pr_pr),...
    sprintf('$T_{total} = %.3f$',c_R_opt_pr_pr'*c_R_opt_pr_pr)},...
    'FontSize',FontSizeVal,'Interpreter','Latex','FitBoxToText','on');
annotation('textbox', [0.2, 0.4, 0.001, 0.001], 'string', ...
    {sprintf('$NB: Ensemble~averaging~is~required$'),
     sprintf('$for~proving~<T_{total}> = 2/3 ~(refer~Code~Package~2)$')},...
    'FontSize',FontSizeVal,'Interpreter','Latex','FitBoxToText','on');
set(gca,'FontSize',18)



%------ Focussing of an evanescent mode focussing via phase conjugation ---
S12_pr_ev=S12(1:init_data.num_modes_prop,(1+init_data.num_modes_prop):end); 
                % Incident modes evanescent, transmitted modes propagating
S21_ev_pr=S21((1+init_data.num_modes_prop):end,1:init_data.num_modes_prop); 

mode_num_to_focus_ev=2;                   % The ith mode to be focussed


I_m_ev=zeros(init_data.num_modes_evanes,1); % Column vector to choose the mode 
                                       %     to focus
I_m_ev(mode_num_to_focus_ev)=1;

c_L_pr_ev=S12_pr_ev*I_m_ev; % outgoing propagating wave on the left 
                            % due the incident single evanescent mode on 
                            % the right 
c1=1/norm(c_L_pr_ev);  % Flux normalization term for phase conjugation to the
                       % right side. 

c_R_opt_ev_ev=S21_ev_pr*(c1.*(conj(S12_pr_ev))*I_m_ev); 
                     % This is the evanescent part of the focussed wave
                     % Phase conjugation to right side of the slab
                     % Here, (c1.*(conj(S12_pr_ev))*I_m_ev) is the incident
                     % flux-normalized phase conjugate wave from the left.  
c_R_opt_pr_ev=S21_pr_pr*(c1.*(conj(S12_pr_ev))*I_m_ev); 
                     % This is the propagating part of the focussed wave
T_background_pr_ev=c_R_opt_pr_ev'*c_R_opt_pr_ev;

figure('Position', [50 50 1200 500],'color','W');
plot(1:init_data.num_modes,[abs(c_R_opt_pr_ev).^2;abs(c_R_opt_ev_ev).^2],'-*')
xlabel('$Mode~number$','Interpreter','Latex')
ylabel('Transmission')
title({'Optimal~transmission through a single disorder : ',...
    'Single~mode~focussing~of~an~\textbf{evanescent}~mode'},...
    'Interpreter','Latex')
annotation('textbox', [0.3, 0.8, 0.001, 0.001], 'string', ...
    {sprintf('$Evanes~mode~num~focussed=%d~(m=%d)$',mode_num_to_focus_ev,...
    init_data.num_modes_prop+mode_num_to_focus_ev),...
    sprintf('$T_{background} = %.3f$',T_background_pr_ev)},...
    'FontSize',FontSizeVal,'Interpreter','Latex','FitBoxToText','on');
annotation('textbox', [0.2, 0.6, 0.001, 0.001], 'string', ...
    {sprintf('$NB: Ensemble~averaging~is~required$'),
     sprintf('$for~proving~<T_{background}> = 2/3 ~(refer~Code~Package~2)$')},...
    'FontSize',FontSizeVal,'Interpreter','Latex','FitBoxToText','on');
set(gca,'FontSize',18)



