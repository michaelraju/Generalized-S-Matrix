size_of_font=20;
figure('Position', [50 50 500 700],'color','W');
colormap(jet);
imagesc([0 init_data.kref*init_data.dz*(size(cascaded_sample,2)-1)],[0 init_data.kref*init_data.dy*(init_data.Ny-1)],(cascaded_sample));  
xlabel('$k_{ref}z$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
axis xy 
axis equal tight
set(gca,'FontSize',size_of_font) 
colorbar
title('$\Delta \epsilon(z,y)$','Interpreter','Latex')
gcf.RendererMode='manual';
%--------------------- Plot the thinner slabs ------------------------------
figure('Position', [0 0 1000 700],'color','W');
subplot(1,4,1)
colormap(bone);
imagesc([0 init_data.kref*init_data.dy*(init_data.Nz-1)],[0 init_data.kref*init_data.dy*(init_data.Ny-1)],eps_profile_array(:,:,1));  
xlabel('$k_{ref}z$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
axis xy 
axis equal tight
set(gca,'FontSize',size_of_font) 
colorbar
title('$\Delta \epsilon(z,y)$','Interpreter','Latex')
gcf.RendererMode='manual';

subplot(1,4,2)
colormap(jet);
imagesc([0 init_data.kref*init_data.dy*(init_data.Nz-1)],[0 init_data.kref*init_data.dy*(init_data.Ny-1)],eps_profile_array(:,:,2));  
xlabel('$k_{ref}z$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
axis xy 
axis equal tight
set(gca,'FontSize',size_of_font) 
colorbar
title('$\Delta \epsilon(z,y)$','Interpreter','Latex')
gcf.RendererMode='manual';

subplot(1,4,3)
colormap(jet);
imagesc([0 init_data.kref*init_data.dy*(init_data.Nz-1)],[0 init_data.kref*init_data.dy*(init_data.Ny-1)],eps_profile_array(:,:,4));  
xlabel('$k_{ref}z$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
axis xy 
axis equal tight
set(gca,'FontSize',size_of_font) 
colorbar
title('$\Delta \epsilon(z,y)$','Interpreter','Latex')
gcf.RendererMode='manual';


subplot(1,4,4)
colormap(jet);
imagesc([0 init_data.kref*init_data.dy*(init_data.Nz-1)],[0 init_data.kref*init_data.dy*(init_data.Ny-1)],eps_profile_array(:,:,end));  
xlabel('$k_{ref}z$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
axis xy 
axis equal tight
set(gca,'FontSize',size_of_font) 
colorbar
title('$\Delta \epsilon(z,y)$','Interpreter','Latex')
gcf.RendererMode='manual';
