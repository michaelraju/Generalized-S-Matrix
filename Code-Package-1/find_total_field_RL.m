function [final_field] = find_total_field_RL(Gij,c_inc,init_data)
%----------------------- initialisation -----------------------------------
kz_flux=init_data.kz_flux;
Y=init_data.Y;
ylin=Y(:,1);
kzdz=init_data.kzdz;
dz=init_data.dz;
Nmat=init_data.Nmat;
W=init_data.W;
Mmat=init_data.Mmat;
kydy=init_data.kydy;
num_modes_prop=init_data.num_modes_prop;
num_modes=init_data.num_modes;
Nz0_R=init_data.Nz0_R;
find_prop_modes_RL=init_data.find_prop_modes_RL;
find_evanes_modes=init_data.find_evanes_modes;

[NY,NZ,~]=size(Gij);
prefact_wave=zeros(NY,NZ);
total_field=zeros(NY,NZ,NY);
final_field=zeros(NY,NZ);

%--------------------------- Propagating modes ----------------------------
for mcount=1:num_modes_prop
prefact_wave=prefact_wave +  ...
    (2*1i*kz_flux(mcount)).*c_inc(mcount).*find_prop_modes_RL(mcount, ...
    kzdz,dz,Nmat,W,Mmat,kydy);
end
%--------------------------- Evanescent modes on the left -----------------
if length(c_inc)> init_data.num_modes_prop
    for mcount=1+num_modes_prop:num_modes
    prefact_wave=prefact_wave + ...
        (-2*kz_flux(mcount)).*c_inc(mcount).*find_evanes_modes(mcount,kzdz, ...
        dz,Nmat,W,Mmat,kydy,Nz0_R);
    end
end  
%--------------------------------------------------------------------------
    for int_count=1:NY
    total_field(:,:,int_count)= ...
        prefact_wave(int_count,Nz0_R).*reshape(Gij(:,:,int_count),NY,NZ);
    end
   
    
    for icount=1:NY
            for jcount=1:NZ
              final_field(icount,jcount)= ...
                  trapz(ylin,flipud(total_field(icount,jcount,:)));
            end
    end
   
end

