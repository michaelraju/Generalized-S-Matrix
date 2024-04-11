FontSizeVal=24;
%------ Focussing of a propagating mode focussing via phase conjugation ---
mode_num_to_focus_pr=10;                   % The ith mode to be focussed
I_m_pr=zeros(init_data.num_modes_prop,1); % Column vector to choose the mode 
                                       %     to focus
I_m_pr(mode_num_to_focus_pr)=1;

clearvars S21_cas S21_big_array S21_prop
c_R_opt_pr_pr=zeros(init_data.num_modes_prop,no_of_thick_slabs);
T_focus_pr_pr=zeros(1,no_of_thick_slabs);
T_background_pr_pr=zeros(1,no_of_thick_slabs);

for scount=1:no_of_thick_slabs
scount    
clearvars S12_cas S21_cas  
load(sprintf('S12_cas_%d.mat',scount));
load(sprintf('S21_cas_%d.mat',scount));

S12_pr_pr=S12_cas(1:init_data.num_modes_prop,1:init_data.num_modes_prop); 
                  % Purely propagating S21 without evanescent modes                  
S21_pr_pr=S21_cas(1:init_data.num_modes_prop,1:init_data.num_modes_prop); 
                  % Purely propagating S21 without evanescent modes 

%--------- Start by exciting the mode to be focussed from the right side
% of the disorder, to be phase conjugated again from the left to focus back
% on the right side. 
c_L_pr_pr=S12_pr_pr*I_m_pr; % outgoing wave on the left due the single mode 
                      % incidence from the right
c1=1/norm(c_L_pr_pr);  % Flux normalization term for phase conjugation to the
                    % right side. 
c_R_opt_pr_pr(:,scount)=S21_pr_pr*(c1.*(conj(S12_pr_pr))*I_m_pr); 
                     %Phase conjugation to right side of the slab
                     % Here, (c1.*(conj(S12_pr_pr))*I_m_pr) is the incident
                     % flux-normalized phase conjugate wave from the left.  

T_focus_pr_pr(scount)=abs(c_R_opt_pr_pr(mode_num_to_focus_pr,scount))^2;
                     % Transmission just due to the focus
T_background_pr_pr(scount)=sum(abs(c_R_opt_pr_pr(1:end ~= mode_num_to_focus_pr,scount)).^2);                    

end

[mean(T_focus_pr_pr) + mean(T_background_pr_pr)]



%------ Focussing of an evanescent mode focussing via phase conjugation ---
mode_num_to_focus_ev=4;                   % The ith mode to be focussed
I_m_ev=zeros(init_data.num_modes_evanes,1); % Column vector to choose the mode 
                                       %     to focus
I_m_ev(mode_num_to_focus_ev)=1;

clearvars S21_cas S21_big_array S21_prop S12_cas S21_pr_pr S12_pr_pr
c_R_opt_ev_ev=zeros(init_data.num_modes_evanes,no_of_thick_slabs);
c_R_opt_pr_ev=zeros(init_data.num_modes_prop,no_of_thick_slabs);
T_background_pr_ev=zeros(1,no_of_thick_slabs);

for scount=1:no_of_thick_slabs
scount    
clearvars S12_cas S21_cas  
load(sprintf('S12_cas_%d.mat',scount));
load(sprintf('S21_cas_%d.mat',scount));

S12_pr_ev=S12_cas(1:init_data.num_modes_prop,(1+init_data.num_modes_prop):end); 
                % Incident modes evanescent, transmitted modes propagating
S21_ev_pr=S21_cas((1+init_data.num_modes_prop):end,1:init_data.num_modes_prop);

S21_pr_pr=S21_cas(1:init_data.num_modes_prop,1:init_data.num_modes_prop); 
                  % Purely propagating S21 without evanescent modes 


c_L_pr_ev=S12_pr_ev*I_m_ev; % outgoing propagating wave on the left 
                            % due the incident single evanescent mode on 
                            % the right 
c1=1/norm(c_L_pr_ev);  % Flux normalization term for phase conjugation to the
                       % right side. 

c_R_opt_ev_ev(:,scount)=S21_ev_pr*(c1.*(conj(S12_pr_ev))*I_m_ev); 
                     % This is the evanescent part of the focussed wave
                     % Phase conjugation to right side of the slab
                     % Here, (c1.*(conj(S12_pr_ev))*I_m_ev) is the incident
                     % flux-normalized phase conjugate wave from the left.  
c_R_opt_pr_ev(:,scount)=S21_pr_pr*(c1.*(conj(S12_pr_ev))*I_m_ev); 
                     % This is the propagating part of the focussed wave
T_background_pr_ev(scount)=c_R_opt_pr_ev(:,scount)'*c_R_opt_pr_ev(:,scount);
end

mean(T_background_pr_ev)


%----------------------- Plotting ----------------------------------------
figure('Position', [50 50 900 600],'color','W');
plot(1:init_data.num_modes_prop,mean(abs(c_R_opt_pr_pr).^2,2),'-o')
xlabel('$Mode~number$','Interpreter','Latex')
ylabel('$Transmission$','Interpreter','Latex')
title({'Optimal~transmission through a single disorder : ',...
    'Single~mode~focusing'},...
    'Interpreter','Latex')
annotation('textbox', [0.2, 0.75, 0.001, 0.001], 'string', ...
    {sprintf('$Prop~mode~num~focused=%d/%d$',mode_num_to_focus_pr,init_data.num_modes_prop),...
    sprintf('$\\langle T_{focus} \\rangle = %.4f$',mean(T_focus_pr_pr)), ...
    sprintf('$\\langle T_{background} \\rangle = %.4f$',mean(T_background_pr_pr)),...
    sprintf('$ \\langle T_{total} \\rangle = \\langle T_{focus} \\rangle+ \\langle T_{background} \\rangle $'),...
    sprintf('$\\qquad \\quad=%.4f \\approx 2/3$',mean(T_focus_pr_pr)+mean(T_background_pr_pr))},...
    'FontSize',FontSizeVal,'Interpreter','Latex','FitBoxToText','on');
    set(gca,'FontSize',FontSizeVal)


hold on
plot(1:init_data.num_modes,[mean(abs(c_R_opt_pr_ev).^2,2);mean(abs(c_R_opt_ev_ev).^2,2)],'-*')
xlabel('$Mode~number$','Interpreter','Latex')
ylabel('Transmission')
title({'Optimal~transmission through a single disorder : ',...
    'Single~mode~focusing'},...
    'Interpreter','Latex')
annotation('textbox', [0.2, 0.3, 0.001, 0.001], 'string', ...
    {sprintf('$Evanes~mode~num~focused=%d~(m=%d)$',mode_num_to_focus_ev,...
    init_data.num_modes_prop+mode_num_to_focus_ev),...
    sprintf('$\\langle T_{background} \\rangle = %.4f \\approx 2/3$',mean(T_background_pr_ev))},...
    'FontSize',FontSizeVal,'Interpreter','Latex','FitBoxToText','on');
legend('$Focusing~on~a~propagating~mode$',...
    '$Focusing~on~an~evanescent~mode$','Interpreter','Latex', 'FontSize',FontSizeVal)
set(gca,'FontSize',FontSizeVal)
set(gca, 'YScale', 'log')

