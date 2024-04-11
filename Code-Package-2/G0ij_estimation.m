function [G] = G0ij_estimation(i0th,j0th,ith,jth,kydy,kzdz,W,num_modes_prop,num_modes,dz)
G=zeros(size(ith));
%------- definition involving propagating modes ---------------------------
for n=1:num_modes_prop
G=G+ (2/W).*sin(kydy(n).*(jth-1)).*sin(kydy(n).*(j0th-1)).*(dz/(2*1i*sin(kzdz(n)))).*(exp(1i.*kzdz(n).*abs(ith-i0th)));
end
%------- definition involving evanescent components -----------------------
for n=(num_modes_prop+1):num_modes
G=G+ (2/W).*sin(kydy(n).*(jth-1)).*sin(kydy(n).*(j0th-1)).*(dz/(-2*sin(abs(kzdz(n))))).*(exp(-abs(kzdz(n)).*abs(ith-i0th)));
end

end

