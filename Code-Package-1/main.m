%---- Code-Package-1 : Generating ensembles of generalized S matrices -----
%---- The code package contains the frequency domain methods for estimating 
%     the generalized scattering matrix with the evanescent wave modes 
%     for a single given disorder.
%---- The code package is also used to visualise various aspects of the
%     wave transport involving shaped waves.
%---- The code is based on the Green's function perturbation method.
%---- Minimum memory requirement is based on the size of the Green's function
%     needed and the perturbation method used, being the direct method 
%     or the indirect method.
%---- The direct matrix inversion method is memory intensive compared 
%     to the iterative method. But could be faster if the sample is
%     thinner.
%---- The Iterative method is slower, but has lower memory requirement and 
%     would be the preferred method for Green's function perturbation for
%     the thicker disorder given in this code package.
%---- One of the two perturbation method can be opted using the variable 
%     'pertubation_method' in the code.
%---- Towards the beginning of the code, free space Green's function needs 
%     to be evaluated from the analytical expression and stored for the source
%     placed in all the particle positions, consuming memory.
%---- The main advantage of such a scheme is that, one could directly
%     estimate the entire Scattering matrix from the perturbed Green's
%     function.
%----- Boundary condition of the scattering problem is the same boundary 
%      condition of the free-space Green's function, being mirror (reflecting)
%      boundary condition (bc) on the transverse boundaries (y axis) and 
%      outgoing bc on the longitudinal boundaries.
clear 
close all

new_run_flag=0;  % Set the flag=1 if a new run is required. 
                 % If flag=1, one may reduce the size of krefW according 
                 % to one's memory availability. 
if new_run_flag==1                 
%----------------- initialise ---------------------------------------------
krefW=400;   % Dimensionless transverse size (along y axis) of the slab
             %               kref*W.    
krefZ=100;    % Dimensionless longitudinal size of the computational domain
             %               kref*Z (along z axis).                           
lambda0=632.8*10^-9;     % Wavelength in nm.
n0=1.9;                  % Reference(Background) medium refractive index.
offset=2; 
  % Offset is defined wrt to the left and right sides of the computational   
  % domain. The Offset value is used to control the dimensionless slab thickness 
  % The dimensionless slab thickness is krefL=krefZ-2*offset,
  % to be initialised later. Larger the offset value for a given krefZ, 
  % thinner the slab.
disorder_corr_flag=1; % 1 for correlated disorder and 0 for uncorrelated.
disorder_spatial_corr_fact=2.0; % A factor used to control disorder
          % spatial correlation length. Larger the factor, larger the 
          % correlation length. Exact value of the dimensionless spatial 
          % correlation length would be estimated later.
          % Usually, disorder_spatial_corr_fact is taken > 1
disorder_perturbation_strength_parameter=0.99;
            % a number between 0 and 1 to adjust the strength of scattering
            % or the degree of the disorder perturbation.
num_modes_evanes=30;  % number of evanescent waves, set by the user.
                      % Larger the number, finer will be the 
                      % spatial resolution, higher the memory 
                      % required for a large sample. 
%pertubation_method='D';  
pertubation_method='Iter';
% 'D' for direct matrix inversion method : Fast, but memory intensive
%      and mainly used for testing purposes involving smaller samples.
% 'Iter' for the iterative method : Slower, but uses lower memory in
%      comparison to the direct matrix inversion method.
block_size_fraction=(1/10); 
% Block size fraction is the fraction of the total no of disorder particle
% grid points taken per step of the iteration process to solve for 
% the perturbed field. For example, if the total no of disorder pertubation 
% grid points is 20,000. then block_size_fraction=(1/5) implies, 
% 4000 particles are solved  at once per the iteration method and the 
% Green's function is updated for every 4000 particles taken together, 
% per iteration step. Hence, 5 steps are needed to update the perturbed field.
init_data = initialisation_frequency_domain(krefW,krefZ,...
    n0,lambda0,num_modes_evanes);  % Initialisation 
init_data.disorder_perturbation_strength_parameter=...
    disorder_perturbation_strength_parameter;
init_data.disorder_corr_flag=disorder_corr_flag;  
init_data.offset=offset;
init_data.krefL=(init_data.Nz-2*init_data.offset-1)*init_data.kref*init_data.dz; 
            % Dimensionless thickness of the 
            % disordered slab along the z axis.
%------------------- define disorder --------------------------------------
init_data.sigma_disorder= ...
  disorder_perturbation_strength_parameter*(init_data.n0^2-1);
                           % dielectric disorder strength
if init_data.disorder_corr_flag==1
init_data.kreflc=disorder_spatial_corr_fact*(init_data.kref*init_data.dz); 
    % kreflc is kref*lc where lc is the spatial correlation length 
    % (init_data.kref*init_data.dz) is the smallest possible value 
    if init_data.krefL<10*init_data.kreflc
      sprintf('Sample thickness is preferred to be atleast 10kref*lc !!!')
        pause;
    end
end
eps_profile = generate_disorder(init_data);
sprintf('No of grid elements = %d',numel(eps_profile))
%-------- Estimate the disorder properties --------------------------------
[kreflc_main_numerical] = ...
    estimate_disorder_spatial_correlation(eps_profile,init_data);
init_data.kreflc_main_numerical=kreflc_main_numerical;
init_data.eps_profile=eps_profile;
init_data.eps_lin_nonzero=find(eps_profile);
init_data.no_of_disorder_perturbations=numel(init_data.eps_lin_nonzero);
init_data.no_of_grid_elements=(init_data.Nz*init_data.Ny);
%------------------- Estimate free space greens function ------------------
if init_data.no_of_disorder_perturbations==0
[G0ij] = evaluate_G0ij(init_data); % Estimating the Green's function 
                                   %    without the disorder.
else 
[G0ik,G0ij] = evaluate_G0ik_G0ij(init_data); % Estimating the Green's 
                                            % function without disorder    
end
G0ij_LR=reshape(G0ij(:,init_data.left_boundary_src), ...
    [init_data.Ny init_data.Nz init_data.Ny]); 
                               % Saving a part of G0ij now to be 
                               %        used for plotting later.
                               % Saved now because in the iterative method, 
                               %         G0ik and G0ij will be overwritten.
%-------------------- Estimate perturbed Green's function -----------------
if init_data.no_of_disorder_perturbations~=0
    if strcmp(pertubation_method,'D') % Choosing the direct method.
    [Gij_R_LR,Gij_L_LR,Gij_LR,Gij_R_RL,Gij_L_RL,Gij_RL]= ...
    greens_function_perturbation_direct(G0ik,G0ij,init_data);
    elseif strcmp(pertubation_method,'Iter')% Choosing the iterative method.
    init_data.block_size_fraction=block_size_fraction;   
    src_ind=[init_data.left_boundary_src init_data.right_boundary_src];
    no_of_iterations=floor(1/init_data.block_size_fraction);
    no_of_particles_per_iteration= ...
        floor(block_size_fraction*init_data.no_of_disorder_perturbations);
    %----------------- Start iteration ------------------------------------
       tic
       for b_count=1:no_of_iterations
        sprintf('Iterative pertubation using scatterer block no %d/%d',...
            b_count,no_of_iterations)
            if b_count~=no_of_iterations    
                particle_index= ...
                1+(b_count-1)*no_of_particles_per_iteration: ...
                b_count*no_of_particles_per_iteration;
            elseif b_count==no_of_iterations
                particle_index= ...
                1+(b_count-1)*no_of_particles_per_iteration: ...
                init_data.no_of_disorder_perturbations;
            end
            % particle_index chooses the block of scatterers taken together 
            % for the perturbing the Green's function
            G0kj=G0ij(init_data.eps_lin_nonzero(particle_index),src_ind);
            % For G0kj, field location k is taken only at the chosen 
            % perturbation particle locations and j the source location at 
            % the boundary.
            G0kk=G0ik(init_data.eps_lin_nonzero(particle_index), ...
                particle_index(end)+1:end); 
            % For G0kk, field location k is taken only at the chosen 
            % perturbation particle locations and the source location k is 
            % taken for the remaining particle locations which are yet to be 
            % brought up for perturbation. 
            %------ Estimate Mpert matrix (Refer the associated paper/thesis)--
            eps_full=(eps_profile(:));    
            k0=init_data.kref/init_data.n0; % Because in Dysons equation, 
                                  % we use the free space k0, not the kref.
            Vkdeltak=...
            -(init_data.dz*init_data.dy*k0^2).*eps_full(init_data.eps_lin_nonzero(particle_index));
            Mpert=(zeros(length(particle_index),length(particle_index))); 
            for scount=1:length(Vkdeltak)
                  for dcount=1:length(Vkdeltak)
                   Mpert(dcount,scount)=...
                       G0ik(init_data.eps_lin_nonzero(particle_index(dcount)),...
                       particle_index(scount))*Vkdeltak(scount);
                  end
            end
            %----------------------- Estimate Gkj -------------------------
            G0kj=(eye(length(Vkdeltak))-Mpert)\G0kj; 
            G0kk=(eye(length(Vkdeltak))-Mpert)\G0kk; 
            clearvars Mpert;
            %--- Estimate Gij (Refer the thesis) --------------------------
            G0ij=G0ij + G0ik(:,particle_index).*repmat(Vkdeltak.',...
                init_data.Nz*init_data.Ny,1)*G0kj;
            temp=G0ik(:,particle_index).*repmat(Vkdeltak.',...
                init_data.Nz*init_data.Ny,1);

            if b_count~=no_of_iterations
                    for r_count=b_count+1:no_of_iterations
                            r_count
                            if r_count~=no_of_iterations    
                                 remaining_particle_index= ...
                                    1+(r_count-1)*no_of_particles_per_iteration: ...
                                    r_count*no_of_particles_per_iteration;
                            elseif r_count==no_of_iterations
                                    remaining_particle_index= ...
                                     1+(r_count-1)*no_of_particles_per_iteration: ...
                                     init_data.no_of_disorder_perturbations;
                            end
                                G0ik(:,remaining_particle_index)= ...
                                    G0ik(:,remaining_particle_index) + ...
                                    temp*G0kk(:,remaining_particle_index-particle_index(end));
                    end
            end

       end
       clearvars G0ik G0kj G0kk temp; % This is a partially overwritten matrix. 
                                 % Good to clear those matrices now.
       %-------- Store the perturbed Green's function matrices ------------
        Gij_L_LR=(G0ij(init_data.left_boundary_field,init_data.left_boundary_src));        
        Gij_R_LR=(G0ij(init_data.right_boundary_field,init_data.left_boundary_src)); 
        Gij_R_RL=(G0ij(init_data.right_boundary_field,init_data.right_boundary_src)); 
        Gij_L_RL=(G0ij(init_data.left_boundary_field,init_data.right_boundary_src)); 
        Gij_LR=reshape(G0ij(:,init_data.left_boundary_src),[init_data.Ny init_data.Nz init_data.Ny]);
        Gij_RL=reshape(G0ij(:,init_data.right_boundary_src),[init_data.Ny init_data.Nz init_data.Ny]);
        % Gij_R_LR : LR means wave incident from left to right, R denotes Green's
        %             function receiver on the right boundary 
        % Gij_L_LR : LR means wave incident from left to right, L denotes Green's
        %             function receiver on the left boundary
        % Gij_LR  : Greens function for the source on the left boundary 
        % Similarly RL means right to left
        sprintf('Time taken iteratively = %f mins',toc/60)
    end
       clearvars G0ij; % Good to clear
      %--------------------------------------------------------------------------    
else 
% If there is no perturbation, then take G=G0, the unperturbed Green's
% functions
Gij_LR=reshape(G0ij(:,init_data.left_boundary_src),[init_data.Ny init_data.Nz init_data.Ny]);
Gij_RL=reshape(G0ij(:,init_data.right_boundary_src),[init_data.Ny init_data.Nz init_data.Ny]);
Gij_L_LR=(G0ij(init_data.left_boundary_field,init_data.left_boundary_src));        
Gij_R_LR=(G0ij(init_data.right_boundary_field,init_data.left_boundary_src)); 
Gij_R_RL=(G0ij(init_data.right_boundary_field,init_data.right_boundary_src)); 
Gij_L_RL=(G0ij(init_data.left_boundary_field,init_data.right_boundary_src));        
end
%-- Estimate the transmission and reflection matrix, S21 and S11 ----------
[S21,S11] = S21S11estimation_generalised(Gij_L_LR.',Gij_R_LR.',init_data);
%-- Estimate Transmission and reflection matrix, S12 and S22 --------------
[S12,S22] = S12S22estimation_generalised(Gij_R_RL.',Gij_L_RL.',init_data);
end

if new_run_flag==0
% If a new run is not needed, load the existing saved data.   
    if isfile('saved-data-Code-Package-1.mat')
     % File exists.
     load('saved-data-Code-Package-1.mat'); 
    else
     % File does not exist and needs to be downloaded from the Zenodo 
     % repository.
     disp('Load the saved-data-Code-Package-1.mat file')
     [selected_data_name,data_path] = uigetfile('*.mat');
         if strcmp(selected_data_name,'saved-data-Code-Package-1.mat')
             load(selected_data_name);
             disp('Loaded saved-data-Code-Package-1.mat file')
         else
             disp('Right .mat not selected !!')
         end
    end   
 
else
disp('Warning : Overwriting the previously saved data...')
pause(5);
save('saved-data-Code-Package-1.mat') % If a new run is needed, save the workspace data 
end


%-------------- Analysis part of the code-package-1  ----------------------
% Validation for the generalised unitarity and reciprocity properties 
generalised_reciprocity_and_unitarity_validation(S11,S12,S21,S22,init_data)          
%------ Script for verifying various flux conservation relations ----------
script_test_flux_conservation
%-- Script for demonstrating the transmission eigenchannel decomposition 
%   involving the propagating part of the S matrix, being S_pr_pr. Script
%   demonstrates wave scattering due to incident shaped waves, exciting 
%   various eigenchannels in the quasi-1D geometry. Incidence from left side
%   or the right side of the disorder is given. 
script_eigenchannels_S21_S12_quasi_1D
%--- If interested in the transmission eigenchannels for the 2D geometry
%    involving a finite sized beam, use the following script.
% script_eigenchannels_S21_finite_sized_beam
%------ Script for testing wave scattering for an unshaped incident wave 
script_test_wave_scattering
%------- Script for visualising single mode focussing ---------------------
script_single_mode_focussing
%---- Additional plots for publication ------------------------------------
script_for_generating_additional_figures
%--------------------------------------------------------------------------