function [S21_lead] = generate_S21_lead(init_data)
%--------------------- Initialisation -------------------------------------
dy=init_data.dy;
dz=init_data.dz;
n0=init_data.n0;
Ny=init_data.Ny;
W=init_data.W;
kydy=init_data.kydy;
kzdz=init_data.kzdz;
kz_flux=init_data.kz_flux(1:init_data.num_modes_prop);
jth=init_data.jth;
num_modes_prop=init_data.num_modes_prop;
kref=init_data.kref;
lambda0=init_data.lambda0;
M_lead_in=ceil(init_data.Ny/init_data.WbyDin);
% z_L is set at 0  Hence N_lead_in=1;
N_lead_in=1;
jth_lead_in=1:M_lead_in;
Chi = @(n,W,jth,kydy) (sqrt(2/W).*sin(kydy(n).*(jth-1)));

if mod(M_lead_in,2)==1
  M_lead_in=M_lead_in+1;  % Make M_lead_in an even number, mainly for the
                          %        simplicity of coding.
  jth_lead_in=1:M_lead_in;
end
sprintf('M_lead_in is %d',M_lead_in)
W_lead_in=(M_lead_in-1)*dy;
%---------------- Find approx no of propagating modes in the lead region --
num_modes_prop_lead_in_approx=ceil(2*W_lead_in*n0/lambda0); % initial guess 
%---------------- Find exact no of propagating modes ----------------------
mode_num_lead_in=1:(num_modes_prop_lead_in_approx+50); % for testing
kydy_lead_in= ((mode_num_lead_in.*pi.*dy)./W_lead_in); % ky * dy 
kzdz_lead_in=  acos( 2-(((kref*dz)^2)/2)-cos(kydy_lead_in) ); %kz * dz    
kz_flux_lead_in= sin(kzdz_lead_in)./dz;
num_modes_prop_lead_in=find(imag(kz_flux_lead_in)==0,1,'last');  
                                 % True estimate of no of propagating modes
if(isreal(kz_flux_lead_in))
   disp('Given number of modes, doesnt span the basis completely')
    pause; 
end
%------------------- Total no of propagating modes ------------------------
num_modes_lead_in=num_modes_prop_lead_in;  % only propagating modes are taken
%------------------- Re-estimate kydy and kzdz ----------------------------
mode_num_lead_in=1:num_modes_lead_in;      
kydy_lead_in= ((mode_num_lead_in.*pi.*dy)./W_lead_in);   %% ky x dy 
kzdz_lead_in=  acos( 2-(((kref*dz)^2)/2)-cos(kydy_lead_in) );  %%  kz x dz    
kz_flux_lead_in= sin(kzdz_lead_in)./dz;
%-------------- Estimate lead diffraction S21 -----------------------------
 S21_lead=zeros(num_modes_prop,num_modes_lead_in); 
                                             % m for input and n for output

 for mcount=1:num_modes_lead_in 
 larger_mode=zeros(Ny,1);
 smaller_mode=(abs(1/sqrt(kz_flux_lead_in(mcount))))*exp(1i*kzdz_lead_in(mcount)*(N_lead_in-1))*Chi(mcount,W_lead_in,jth_lead_in,kydy_lead_in).'; 
 larger_mode(1+floor(Ny/2):floor(Ny/2)+M_lead_in/2)=smaller_mode(1+M_lead_in/2:M_lead_in);
 larger_mode(1+floor(Ny/2)-M_lead_in/2:floor(Ny/2))=smaller_mode(1:M_lead_in/2);
for ncount=1:num_modes_prop
  chi_m=Chi(ncount,W,jth,kydy).';
  S21_lead(ncount,mcount)=abs(sqrt(kz_flux(ncount)))*sum(chi_m.*larger_mode.*dy);
end
end
 
%---------------- Test for flux convervation ------------------------------
S21_lead_prop=S21_lead(1:num_modes_prop,1:num_modes_prop_lead_in); 
trans_diff=trace(S21_lead_prop'*S21_lead_prop)./num_modes_prop_lead_in;
sprintf('S21_lead transmission=%f',trans_diff)

%figure
for mcount=1:num_modes_lead_in %% input lead
  mcount  
  smaller_mode=(abs(1/sqrt(kz_flux_lead_in(mcount))))*exp(1i*kzdz_lead_in(mcount)*(N_lead_in-1))*Chi(mcount,W_lead_in,jth_lead_in,kydy_lead_in).'; 
  larger_mode(1+floor(Ny/2):floor(Ny/2)+M_lead_in/2)=smaller_mode(1+M_lead_in/2:M_lead_in);
  larger_mode(1+floor(Ny/2)-M_lead_in/2:floor(Ny/2))=smaller_mode(1:M_lead_in/2);
  chi_m=zeros(length(jth),1);
  for ncount=1:num_modes_prop
  chi_m=chi_m+S21_lead(ncount,mcount).*(abs(1/sqrt(kz_flux(ncount)))).*exp(1i*kzdz(ncount)*(N_lead_in-1)).*Chi(ncount,W,jth,kydy).';
  end 
%-------------- Plot S21 diffraction --------------------------------------
%   plot(abs(chi_m),'-')
%   hold on
%   plot(abs(larger_mode),'*')
%   title(sprintf('Mode no %d',mcount));
%   legend('Fit','Original')
%   hold off
%   drawnow
%   pause(0.2)
end
%--------------------------------------------------------------------------
figure('position',[100 100 1600 600])
subplot(1,3,1)
colormap jet
imagesc(real(S21_lead))
axis equal tight
colorbar
xlabel('$incident~mode~no$','Interpreter','Latex')
ylabel('$transmitted~mode~no$','Interpreter','Latex')
title('$real(S^{lead}_{21})$','Interpreter','Latex')
set(gca,'FontSize',18) 

subplot(1,3,2)
colormap jet
imagesc(imag(S21_lead))
axis equal tight
xlabel('$incident~mode~no$','Interpreter','Latex')
ylabel('$transmitted~mode~no$','Interpreter','Latex')
title('$imag(S^{lead}_{21})$','Interpreter','Latex')
set(gca,'FontSize',18) 
colorbar

subplot(1,3,3)
colormap jet
imagesc(abs(S21_lead'*S21_lead))
axis equal tight
title('${S^{lead}_{21}}^\dagger S^{lead}_{21}$','Interpreter','Latex')
set(gca,'FontSize',18) 
colorbar
%----------------------- Annotation ---------------------------------------
side_description={'$No~of~finite$',...
                  '$sized~propagating$', ...
    sprintf('$incident~modes$ = %d',num_modes_prop_lead_in),...
    sprintf('$W/D_{in}=%.3f$',init_data.WbyDin), ...
    sprintf('$S_{21}^{lead} transmission$'), ...
    sprintf('$=%.4f$',trans_diff)};
annotation('textbox', [0.005, 0.95, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on');   
%--------------------------------------------------------------------------
end

