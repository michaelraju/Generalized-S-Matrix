function [Gij_R_LR,Gij_L_LR,Gij_LR,Gij_R_RL,Gij_L_RL,Gij_RL] = greens_function_perturbation_direct(G0ik,G0ij,init_data)
%-------------------- Initialisation --------------------------------------
eps_profile=init_data.eps_profile;
Ny=init_data.Ny;
Nz=init_data.Nz;
eps_lin_nonzero=init_data.eps_lin_nonzero;
kref=init_data.kref;
n0=init_data.n0;
dy=init_data.dy;
dz=init_data.dz;
left_boundary_src=init_data.left_boundary_src;
right_boundary_src=init_data.right_boundary_src;
left_boundary_field=init_data.left_boundary_field;
right_boundary_field=init_data.right_boundary_field;
src_ind=[left_boundary_src right_boundary_src];
%-------------- Estimate G0kj only at particle positions ------------------
G0kj=(G0ij(eps_lin_nonzero,src_ind)); %% 
%----------------- Estimate Mpert matrix (Refer thesis) -------------------
eps_full=(eps_profile(:));    
k0=kref/n0;   % Because in Dysons equation,
              %      we use the free space k0, not the kref.
Vkdeltak=-(dz*dy*k0^2).*eps_full(eps_lin_nonzero);
Mpert=(zeros(length(eps_lin_nonzero),length(eps_lin_nonzero))); %% 
tic
sprintf('Evaluating the M matrix in the Dysons equation... ')
for scount=1:length(Vkdeltak)
    for dcount=1:length(Vkdeltak)
    Mpert(dcount,scount)=G0ik(eps_lin_nonzero(dcount),scount)*Vkdeltak(scount);
    end
end
%----------------------- Estimate Gkj ------------------------------------
sprintf('Evaluating Gkj_LR in the Dysons equation...')
Gkj=(eye(length(Vkdeltak))-Mpert)\G0kj; 
clearvars Mpert;
%--- Estimate Gij (Refer the thesis) -------------------------------------
Gij=G0ij + G0ik.*repmat(Vkdeltak.',Nz*Ny,1)*Gkj;
%----------------- Return the perturbed Green's function matrices ---------
Gij_L_LR=(Gij(left_boundary_field,left_boundary_src));        
Gij_R_LR=(Gij(right_boundary_field,left_boundary_src)); 
Gij_R_RL=(Gij(right_boundary_field,right_boundary_src)); 
Gij_L_RL=(Gij(left_boundary_field,right_boundary_src)); 
%--------------------------------------------------------------------------
Gij_LR=reshape(Gij(:,init_data.left_boundary_src),[init_data.Ny init_data.Nz init_data.Ny]);
Gij_RL=reshape(Gij(:,init_data.right_boundary_src),[init_data.Ny init_data.Nz init_data.Ny]);
sprintf('Time taken iteratively = %f mins',toc/60)
%--------------------------------------------------------------------------
end

