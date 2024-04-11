%------------------- define disorder --------------------------------------
init_data.disorder_corr_flag=1;  % set 0 for spatially uncorrelated disorder 

if init_data.disorder_corr_flag==1
init_data.kreflc=kreflc_factor*(init_data.kref*init_data.dz);
    % kreflc is kref*lc where lc is the spatial correlation length 
    % (init_data.kref*init_data.dz) is the smallest possible value 
    if init_data.krefL<5*init_data.kreflc
        sprintf('Sample thickness should be atleast 5kref*lc !!!')
        pause;
    end
end

%---------------- generate S matrices -------------------------------------
S21_array=zeros(init_data.num_modes,init_data.num_modes,no_of_samples);
S11_array=zeros(init_data.num_modes,init_data.num_modes,no_of_samples);
S12_array=zeros(init_data.num_modes,init_data.num_modes,no_of_samples);
S22_array=zeros(init_data.num_modes,init_data.num_modes,no_of_samples);
eps_profile_array=zeros(init_data.Ny,init_data.Nz,no_of_samples);
init_data.kreflc_main_numerical_array=zeros(1,no_of_samples);
init_data.disorder_var_array=zeros(1,no_of_samples);

for ens_count=1:no_of_samples
sprintf('Evaluating sample no %d',ens_count)
%-------- estimate the disorder and its properties ------------------------
eps_profile = generate_disorder(init_data);
sprintf('No of grid elements = %d',numel(eps_profile))
%-------- Estimate the disorder properties --------------------------------
[kreflc_main_numerical,var_disorder] = ...
    estimate_disorder_spatial_correlation(eps_profile,init_data);
init_data.kreflc_main_numerical_array(ens_count)=kreflc_main_numerical;
init_data.disorder_var_array(ens_count)=var_disorder;
sprintf('Dimensionless spatial correlation length kref*lc= %f',init_data.kreflc_main_numerical_array(ens_count))

init_data.eps_profile=eps_profile;
init_data.eps_lin_nonzero=find(eps_profile);
init_data.no_of_disorder_perturbations=numel(init_data.eps_lin_nonzero);
init_data.no_of_grid_elements=(init_data.Nz*init_data.Ny);
%-------------------- Estimate perturbed Green's function -----------------
[Gij_R_LR,Gij_L_LR,Gij_LR,Gij_R_RL,Gij_L_RL,Gij_RL]= ...
    greens_function_perturbation_direct(G0ik,G0ij,init_data);% G1ij_R_LR : LR means wave incident from left to right, R denotes Green's
%             function receiver on the right boundary 
% Gij_L_LR : LR means wave incident from left to right, L denotes Green's
%             function receiver on the left boundary
% Gij_LR  : Greens function for the source on the left boundary 
% Similarly RL means right to left
%-- Estimate the transmission and reflection matrix, S21 and S11 ----------
[S21,S11] = S21S11estimation_generalised(Gij_L_LR.',Gij_R_LR.',init_data);
%-- Estimate Transmission and reflection matrix, S12 and S22 --------------
[S12,S22] = S12S22estimation_generalised(Gij_R_RL.',Gij_L_RL.',init_data);

S21_array(:,:,ens_count)=S21;
S11_array(:,:,ens_count)=S11;
S12_array(:,:,ens_count)=S12;
S22_array(:,:,ens_count)=S22;
eps_profile_array(:,:,ens_count)=eps_profile;
pause(1);
close all;
end
%----------------------------- Save data ----------------------------------
% save('S21_array.mat','S21_array');
% save('S11_array.mat','S11_array');
% save('S12_array.mat','S12_array');
% save('S22_array.mat','S22_array');
% save('eps_profile_array.mat','eps_profile_array');
init_data.no_of_samples=no_of_samples;
% save('init_data.mat','init_data');
