function [mask_prop_prop,mask_evan_evan,mask_evan_prop,mask_prop_evan] =...
    generate_selection_mask(num_modes,num_modes_prop)
FontSizeVal=18;
num_modes_evanes=num_modes-num_modes_prop;
%----------------- initialisation -----------------------------------------
mask_prop_prop=zeros(2*num_modes,2*num_modes);
mask_evan_evan=zeros(2*num_modes,2*num_modes);
mask_evan_prop=zeros(2*num_modes,2*num_modes);
%--------------------------------------------------------------------------
mask_prop_prop([1:num_modes_prop num_modes+1:num_modes+num_modes_prop],[1:num_modes_prop num_modes+1:num_modes+num_modes_prop])=1;
mask_evan_evan(num_modes_prop+1:num_modes,num_modes_prop+1:num_modes)=1;
mask_evan_evan(2*num_modes-num_modes_evanes+1:2*num_modes,2*num_modes-num_modes_evanes+1:2*num_modes)=1;
mask_evan_evan(num_modes_prop+1:num_modes,num_modes_prop+1:num_modes)=1;
mask_evan_evan(2*num_modes-num_modes_evanes+1:2*num_modes,num_modes_prop+1:num_modes)=1;
mask_evan_evan(num_modes_prop+1:num_modes,2*num_modes-num_modes_evanes+1:2*num_modes)=1;

mask_evan_prop(num_modes_prop+1:num_modes,1:num_modes_prop)=1;
mask_evan_prop(2*num_modes-num_modes_evanes+1:2*num_modes,1:num_modes_prop)=1;
mask_evan_prop(num_modes_prop+1:num_modes,(num_modes+1):(num_modes+num_modes_prop))=1;
mask_evan_prop(2*num_modes-num_modes_evanes+1:2*num_modes,(num_modes+1):(num_modes+num_modes_prop))=1;

mask_prop_evan=-(mask_prop_prop+mask_evan_evan+mask_evan_prop-1);
%--------------------------------------------------------------------------
figure('Position', [50 50 800 800],'color','W');
subplot(2,2,1)
spy(mask_prop_prop)
xlabel('');
title('$Mask(pr,pr)$','Interpreter','Latex')
xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes+1) (num_modes+floor(num_modes_prop/4)) (num_modes+floor(num_modes_prop/2)) (num_modes+floor(3*num_modes_prop/4)) (num_modes+num_modes_prop) 2*num_modes]);
ticks= [1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop            1             floor(num_modes_prop/4)             floor(num_modes_prop/2)             floor(3*num_modes_prop/4)             num_modes_prop   num_modes_evanes];
yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes+1) (num_modes+floor(num_modes_prop/4)) (num_modes+floor(num_modes_prop/2)) (num_modes+floor(3*num_modes_prop/4)) (num_modes+num_modes_prop) 2*num_modes]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
set(gca,'FontSize',FontSizeVal)

subplot(2,2,2)
spy(mask_evan_evan)
xlabel('');
title('$Mask(ev,ev)$','Interpreter','Latex')
xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes+1) (num_modes+floor(num_modes_prop/4)) (num_modes+floor(num_modes_prop/2)) (num_modes+floor(3*num_modes_prop/4)) (num_modes+num_modes_prop) 2*num_modes]);
ticks= [1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop            1             floor(num_modes_prop/4)             floor(num_modes_prop/2)             floor(3*num_modes_prop/4)             num_modes_prop   num_modes_evanes];
yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes+1) (num_modes+floor(num_modes_prop/4)) (num_modes+floor(num_modes_prop/2)) (num_modes+floor(3*num_modes_prop/4)) (num_modes+num_modes_prop) 2*num_modes]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
set(gca,'FontSize',FontSizeVal)


subplot(2,2,3)
spy(mask_evan_prop)
xlabel('');
title('$Mask(ev,pr)$','Interpreter','Latex')
xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes+1) (num_modes+floor(num_modes_prop/4)) (num_modes+floor(num_modes_prop/2)) (num_modes+floor(3*num_modes_prop/4)) (num_modes+num_modes_prop) 2*num_modes]);
ticks= [1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop            1             floor(num_modes_prop/4)             floor(num_modes_prop/2)             floor(3*num_modes_prop/4)             num_modes_prop   num_modes_evanes];
yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes+1) (num_modes+floor(num_modes_prop/4)) (num_modes+floor(num_modes_prop/2)) (num_modes+floor(3*num_modes_prop/4)) (num_modes+num_modes_prop) 2*num_modes]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
set(gca,'FontSize',FontSizeVal)


subplot(2,2,4)
spy(mask_prop_evan)
xlabel('');
title('$Mask(pr,ev)$','Interpreter','Latex')
xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes+1) (num_modes+floor(num_modes_prop/4)) (num_modes+floor(num_modes_prop/2)) (num_modes+floor(3*num_modes_prop/4)) (num_modes+num_modes_prop) 2*num_modes]);
ticks= [1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop            1             floor(num_modes_prop/4)             floor(num_modes_prop/2)             floor(3*num_modes_prop/4)             num_modes_prop   num_modes_evanes];
yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes+1) (num_modes+floor(num_modes_prop/4)) (num_modes+floor(num_modes_prop/2)) (num_modes+floor(3*num_modes_prop/4)) (num_modes+num_modes_prop) 2*num_modes]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
set(gca,'FontSize',FontSizeVal)

%--------------------------- Test -----------------------------------------
if isempty(find((mask_prop_prop+mask_evan_evan+mask_evan_prop+mask_prop_evan)-ones(2*num_modes,2*num_modes)))==0
    disp('Error !!!')
end
end

