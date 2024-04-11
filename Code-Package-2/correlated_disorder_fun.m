function [corr_disorder_main_uniform]=correlated_disorder_fun(init_data)
%------------------------ Initialisation ----------------------------------
kref=init_data.kref;
sigma_disorder=init_data.sigma_disorder;
Ny=init_data.Ny;
Nz=init_data.Nz;
dy=init_data.dy;
dz=init_data.dz;
kreflc=init_data.kreflc;
rng('shuffle');    %  random number seed based on time 
eps_profile=sigma_disorder.*(2.*rand(Ny,Nz)-1); 
%---------- generate the sample -------------------------------------------
N_lc_main_analyt=ceil(kreflc/(kref*dz)); % Correlation length in grid units
normcdf_fun= @(x) 1/2*erfc(-x/sqrt(2)); % CDF for a normalised gaussian
cdf_fun = @(x,W_uniform) ((0.5.*W_uniform).*(2.*normcdf_fun(x)-1));   
                % A function for converting a normal distribution to a 
                % uniform distribution.
W_white=max(eps_profile(:))-min(eps_profile(:));
corr_disorder_main_normalised = generate_correlated_disorder(eps_profile,kref,dy,dz,N_lc_main_analyt);
corr_disorder_main_uniform=cdf_fun(corr_disorder_main_normalised,W_white);
end

function [corr_disorder_normalised] = generate_correlated_disorder(eps_profile,kref,dy,dz,Nlc)
[Ny,Nz]=size(eps_profile);
%------------------- Step 2 : Define Gaussian power spectrum --------------
max_y=kref*dy*(Ny-1);  % max distance taken along +y axis
max_z=kref*dz*(Nz-1);  % max distance taken along +z axis

ds=(kref*dz);      % space discretisation step 
lc=Nlc*(kref*dz);
y = 0:ds:max_y;    % space vector 
z=  0:ds:max_z;    % space vector 
[z_mesh,y_mesh]=meshgrid(z,y);
wy = lc/sqrt(2);   % so that wy^2 + wz^2 = lc^2
wz = lc/sqrt(2);   
A = 1;             % parameters for the gaussian function
ycen=max_y/2;
zcen=max_z/2;
gauss_corr = A*exp(((-1).*((y_mesh-ycen)./wy).^2) + ...
            ((-1).*((z_mesh-zcen)./wz).^2) ); % Gaussian function
%------ step 3 : Take fourier transform of uncorrelated disorder ----------
disorder_fft=fftshift(fft2(eps_profile));
fft_corr=fftshift(fft2(gauss_corr));
%------- Step 4 : Take the product of disorder_fft and powerspectrum ------
conv_prod=disorder_fft.*fft_corr;
%-------- Step 5 : Estimate Inverse FFT -----------------------------------
corr_disorder = ifft2(ifftshift(conv_prod));  % FFT transform
%-------- Step 6 : Normalise the correlated disorder for unit 
%         standard deviation                         
corr_disorder_normalised=(corr_disorder)./std(corr_disorder(:));
sprintf('Standard deviation of the corrrelated sample is %f', ...
    std(corr_disorder_normalised(:)))
end

