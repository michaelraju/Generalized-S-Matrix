% Plot options are 'real_part','imaginary_part','magnitude'
%plot_type='real_part';  
%plot_type='imaginary_part';  
plot_type='magnitude';  
c_inc=zeros(init_data.num_modes,1);
inc_mode_no=3;          % Incident propagating eigenmode no 
c_inc(inc_mode_no)=1;   % Coefficient of the incident wave
% The incident mode could also be evanescent as the following 
% c_inc(init_data.num_modes_prop+1)=1; 
% The incident mode could be a random propagating wave as the following
% c_inc(1:init_data.num_modes_prop)=(2.*rand(init_data.num_modes_prop,1)-1)+ ...
%                                  1i.*(2.*rand(init_data.num_modes_prop,1)-1);
% c_inc=c_inc./norm(c_inc); % flux normalization

if (~isempty(find(c_inc(1:init_data.num_modes_prop),1)))
c_inc=c_inc./norm(c_inc(1:init_data.num_modes_prop));   
                   % Ensure that the flux of the propagating part of the 
                   %            incident wave has to be normalised to unity
end
field_visualisation_and_comparison_LR(c_inc,Gij_LR,S21,S11,init_data,plot_type)
side_description={'$Wave~scattering$','$due~to~an$','$incident~wave$'};
annotation('textbox', [0.005, 0.95, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on');   



%--- Testing an incident wavefront : propagation from right to left
c_inc=zeros(init_data.num_modes,1);
c_inc(inc_mode_no)=3;
% The incident mode could also be evanescent as the following 
% c_inc(init_data.num_modes_prop+1)=1; 
% The incident mode could be a random propagating wave as the following
% c_inc(1:init_data.num_modes_prop)=(2.*rand(init_data.num_modes_prop,1)-1)+ ...
%                                  1i.*(2.*rand(init_data.num_modes_prop,1)-1);
% c_inc=c_inc./norm(c_inc); % flux normalization

if (~isempty(find(c_inc(1:init_data.num_modes_prop),1)))
c_inc=c_inc./norm(c_inc(1:init_data.num_modes_prop));   
                   % Flux of the propagating part of the 
                   %            incident wave has to be normalised to unity
end
field_visualisation_and_comparison_RL(c_inc,Gij_RL,S12,S22,init_data,plot_type)
side_description={'$Wave~scattering$','$due~to~an$','$incident~wave$'};
annotation('textbox', [0.005, 0.95, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on');  