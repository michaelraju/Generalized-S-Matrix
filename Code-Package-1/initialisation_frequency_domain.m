% Frequency domain method for estimating the generalised scattering matrix 
%      with evanescent waves for a given disorder.
% Main purpose of this code is to visualise various aspects of the
%      framework, discussed in chapter 1 of the thesis.
%

function [init_data] = initialisation_frequency_domain(krefW,krefZ,n0,lambda0,num_modes_evanes)
%-------------------- initialisation---------------------------------------
kref=2*pi*n0/lambda0;      % Reference(Background) medium wave vector.

W=krefW/kref;              % Transverse size of the slab, W
Zmax=krefZ/kref;            % Longitudinal size of the domain, Zmax

num_modes_prop_approx=ceil(2*W*n0/lambda0);  % Approximate no of prop modes 
                                            %     as an initial guess.

% Choose grid resolution 
if num_modes_evanes==0
dy=1/kref;       % grid resolution dy, such that kref*dy=1.
else
kymax=(num_modes_prop_approx+num_modes_evanes)*pi/W; % Largest transverse 
                                                     % spatial frequency.
fact=kref/kymax;% Define the grid resolution factor for having kymax*dy=1.
                %   so that grid dy reduces for incorporating 
                %    evanescent modes.      
dy=fact/kref;   % Grid resolution dy.
end

dz=dy;          % Always keep dy=dz for the validity of the analytical 
                % formulations used.
Nz=floor(Zmax/dz) + 1;  % grid size along z.
Ny=floor(W/dy) + 1;     % grid size along y.
 
W=(Ny-1)*dy;     % snap the grid to fit integer dy.
Zmax=(Nz-1)*dz;  % snap the grid to fit integer dz.

Nz0_L=1; % Position index where evanescent waves begin (left side of the slab)
Nz0_R=Nz; % Position index where evanescent waves begin (Right side of the slab)

% Generate a mesh
jth=1:Ny;
[Nmat,Mmat]=meshgrid(1:Nz,1:Ny);    
Z=(Nmat-1).*dz;
Y=(Mmat-1).*dy;

jind=[1:Ny Ny*Nz-Ny+1:Ny*Nz]; % Linear source indices for the boundary region.
%---------------- find exact no of propagating modes ----------------------
% Testing
mode_num=1:(num_modes_prop_approx+50);       % 50 just for testing purposes
kydy= ((mode_num.*pi.*dy)./W);               % ky*dy 
kzdz= acos( 2-(((kref*dz)^2)/2)-cos(kydy) ); % kz*dz    
kz_flux= sin(kzdz)./dz; % kz_flux is re-estimated later
num_modes_prop=find(imag(kz_flux)==0,1,'last');   
                                 % True estimate of no of propagating modes

% True estimates by reestimating 
num_modes=num_modes_evanes+num_modes_prop;  % total number of modes
%------------------- re-estimate kydy and kzdz ----------------------------
mode_num=1:num_modes;         
kydy= ((mode_num.*pi.*dy)./W);  % True estimate, ky*dy 
kzdz= acos( 2-(((kref*dz)^2)/2)-cos(kydy) );  % True estimate,  kz*dz    
kz_flux= sin(abs(kzdz))./dz;    % kz_flux is re-estimated

% Tests
if (dy~=dz)
 sprintf('Error : dy should be = dz for the analytic formulation used in the code to be valid')
 pause;
end

if(isreal(kzdz))
   disp('Given number of modes, doesnt completely span the basis (including defined no of evanescent modes) ')
    pause; 
end

%----------------------- Define  basis functions -------------------------
 Chi = @(n,W,jth,kydy) (sqrt(2/W).*sin(kydy(n).*(jth-1)));
%--------------------- Propagating part ----------------------------------
phi_LR_prop = @(n,kzdz,dz,Nmat) sqrt(dz./sin(kzdz(n))).*(exp(+1i.*kzdz(n).*(Nmat-1)));
phi_RL_prop = @(n,kzdz,dz,Nmat) sqrt(dz./sin(kzdz(n))).*(exp(-1i.*kzdz(n).*(Nmat-1)));
find_prop_modes_LR =@(mcount,kzdz,dz,Nmat,W,Mmat,kydy) (phi_LR_prop(mcount,kzdz,dz,Nmat).*Chi(mcount,W,Mmat,kydy));
find_prop_modes_RL =@(mcount,kzdz,dz,Nmat,W,Mmat,kydy) (phi_RL_prop(mcount,kzdz,dz,Nmat).*Chi(mcount,W,Mmat,kydy));
%--------------------- Evanescent  part is coordinate dependent ----------
phi_evanes = @(n,kzdz,dz,Nmat,Nz0)  sqrt(dz./sin(abs(kzdz(n)))).*exp(-abs(kzdz(n)).*abs((Nmat-1)-(Nz0-1))) ;  % defined wrt Nz_0
find_evanes_modes=@(mcount,kzdz,dz,Nmat,W,Mmat,kydy,Nz0) (phi_evanes(mcount,kzdz,dz,Nmat,Nz0).*Chi(mcount,W,Mmat,kydy));
% Test orthogonality 
 sum(Chi(num_modes,W,jth,kydy).*Chi(num_modes,W,jth,kydy).*dz);   %%% orthogonality property
 sum(Chi(num_modes,W,jth,kydy).*Chi(num_modes-1,W,jth,kydy).*dz); %%% orthogonality property
%---------------------- Save the structure containing intialisations ------
init_data.krefW=(Ny-1)*kref*dy;
init_data.krefZ=(Nz-1)*kref*dz;
init_data.lambda0=lambda0;
init_data.n0=n0;
init_data.kref=kref;
init_data.W=W;
init_data.Zmax=Zmax;
init_data.num_modes_prop=num_modes_prop;
init_data.num_modes_evanes=num_modes_evanes;
init_data.num_modes=num_modes;
init_data.dy=dy;
init_data.dz=dz;
init_data.kydy=kydy;
init_data.kzdz=kzdz;
init_data.kz_flux=kz_flux;
init_data.kz_flux=kz_flux;
init_data.Z=Z;
init_data.Y=Y;
init_data.jth=jth;
init_data.Nmat=Nmat;
init_data.Mmat=Mmat;
init_data.Nz=Nz;
init_data.Ny=Ny;

init_data.Chi=Chi;
init_data.phi_LR_prop=phi_LR_prop; 
init_data.phi_RL_prop=phi_RL_prop;
init_data.find_prop_modes_LR=find_prop_modes_LR;
init_data.find_prop_modes_RL=find_prop_modes_RL;

init_data.phi_evanes=phi_evanes;
init_data.find_evanes_modes=find_evanes_modes;

init_data.jind=jind; 
init_data.left_boundary_src=1:init_data.Ny;
init_data.right_boundary_src=init_data.Ny+1:2*init_data.Ny;
init_data.left_boundary_field=1:init_data.Ny;
init_data.right_boundary_field=init_data.Ny*init_data.Nz- ...
    init_data.Ny+1:init_data.Ny*init_data.Nz;
init_data.Nz0_L=Nz0_L;  % z index for the left boundary
init_data.Nz0_R=Nz0_R;  % z index for the right boundary
init_data.num_modes_evanes=num_modes_evanes;
%--------------------------------------------------------------------------
end

