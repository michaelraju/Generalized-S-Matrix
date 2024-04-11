function [G0ik,G0ij] = evaluate_G0ik_G0ij(init_data)
%--------------------- initialisation -------------------------------------
Ny=init_data.Ny;
Nz=init_data.Nz;
Nmat=init_data.Nmat;
Mmat=init_data.Mmat;
kydy=init_data.kydy;
kzdz=init_data.kzdz;
W=init_data.W;
num_modes_prop=init_data.num_modes_prop;
num_modes=init_data.num_modes;
dz=init_data.dz;
eps_lin_nonzero=init_data.eps_lin_nonzero;
no_of_disorder_perturbations=init_data.no_of_disorder_perturbations;

jind=init_data.jind; 
%-------------------------- Evaluate G0ik ---------------------------------
G0ik=(zeros(Ny*Nz,no_of_disorder_perturbations)); 
                        % source position indexed using k, where k denotes
                        %    the location at which perturbation is supposed
                        %    to be brought in.
                        % Field evaluation position indexed i.
                        % source taken all over the scatterer position.
                        
G0ij=(zeros(Ny*Nz,2*Ny)); % source position indexed using j.
                          % Field evaluation position indexed i.
                          % source taken over the left and 
                          %       the right boundaries.                        
                        
tic
sprintf('Evaluating free space Greens function G0ik ...')
for kcount=1:no_of_disorder_perturbations    
sprintf('Evaluating G0ik where the source is at disorder location no %d/%d', ...
    kcount,no_of_disorder_perturbations)    
[nj_y,nj_z]=ind2sub([Ny Nz],eps_lin_nonzero(kcount));
[G] =  G0ij_estimation(nj_z,nj_y,Nmat,Mmat,kydy,kzdz,W,...
    num_modes_prop,num_modes,dz);
G0ik(:,kcount)=G(:);  % source index taken as the last index
                      %          variable within the bracket. 
end
sprintf('Time taken to evaluate G0ik : %f mins',toc/60)
%------------------------ Evaluate G0ij -----------------------------------
tic
sprintf('Evaluating free space Greens function G0ij ...')
for jcount=1:2*Ny    
sprintf('Evaluating G0ij where the source is at boundary location no %d/%d', ...
    jcount,2*Ny)    
[nj_y,nj_z]=ind2sub([Ny Nz],jind(jcount));
[G] =  G0ij_estimation(nj_z,nj_y,Nmat,Mmat,kydy,kzdz,W,num_modes_prop, ...
    num_modes,dz);
G0ij(:,jcount)=G(:);  % source index taken as the last index variable 
                      %     within the bracket. 
end
sprintf('Time taken to evaluate G0ij : %f mins',toc/60)
%--------------------------------------------------------------------------
end

