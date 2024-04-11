close all
%------------------ Cascading parameters ----------------------------------
krefW=909;
no_of_samples=25;
num_modes_evanes=30;
lambda0=632.8*10^-9;     % Wavelength in nm.
n0=1.6;                  % Reference(Background) medium refractive index.
disorder_perturbation_strength_parameter=0.65; % kept less than 1
krefZ=15; 
offset=1; 
kreflc_factor=[1];
%--------------------------------------------------------------------------
[init_data] = initialisation_frequency_domain(krefW,krefZ,n0,lambda0,num_modes_evanes);
init_data.disorder_perturbation_strength_parameter=...
    disorder_perturbation_strength_parameter;
init_data.offset=offset;
init_data.krefL=(init_data.Nz-2*init_data.offset-1)*init_data.kref*init_data.dz; 
            % Dimensionless thickness of the 
            % disordered slab along the z axis.
%------ generate a sample disorder (uncorrelated) -----------
init_data.sigma_disorder= ...
  disorder_perturbation_strength_parameter*(init_data.n0^2-1);
init_data.disorder_corr_flag=0; % 1 for correlated disorder and 0 for uncorrelated.
eps_profile = generate_disorder(init_data);
init_data.eps_profile=eps_profile;
init_data.eps_lin_nonzero=find(eps_profile);
init_data.no_of_disorder_perturbations=numel(init_data.eps_lin_nonzero);
init_data.no_of_grid_elements=(init_data.Nz*init_data.Ny);
%------------- Estimate unperturbed Green's function ----------------------
[G0ik,G0ij] = evaluate_G0ik_G0ij(init_data); % Estimating the Green's 
                                            % function without disorder    
