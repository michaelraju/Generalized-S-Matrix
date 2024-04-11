function generalised_reciprocity_and_unitarity_validation(S11,S12,S21,S22,init_data)
Smatrix=[S11 S12; S21 S22]; % Combining to form an S matrix
FontSizeVal=18;
% First plot the S matrix obtained from the generalised Fisher-Lee
%       relations
num_modes_prop=init_data.num_modes_prop;
num_modes=init_data.num_modes;
num_modes_evanes=init_data.num_modes_evanes;

figure('Position', [50 50 1500 500],'color','W');
subplot(1,2,1)
colormap jet
imagesc((abs(Smatrix).^2))
axis equal tight
colorbar
xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes+1) (num_modes+floor(num_modes_prop/4)) (num_modes+floor(num_modes_prop/2)) (num_modes+floor(3*num_modes_prop/4)) (num_modes+num_modes_prop) 2*num_modes]);
ticks= [1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop            1             floor(num_modes_prop/4)             floor(num_modes_prop/2)             floor(3*num_modes_prop/4)             num_modes_prop   num_modes_evanes];
yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes+1) (num_modes+floor(num_modes_prop/4)) (num_modes+floor(num_modes_prop/2)) (num_modes+floor(3*num_modes_prop/4)) (num_modes+num_modes_prop) 2*num_modes]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
title('$|S|^2$','Interpreter','Latex')
set(gca,'FontSize',FontSizeVal)
caxis([min(min((abs(Smatrix).^2))) max(max((abs(Smatrix).^2)))])
set(gca,'ColorScale','log')

subplot(1,2,2)
colormap jet
imagesc(angle(Smatrix))
axis equal tight
caxis([-pi pi])
colorbar
xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes+1) (num_modes+floor(num_modes_prop/4)) (num_modes+floor(num_modes_prop/2)) (num_modes+floor(3*num_modes_prop/4)) (num_modes+num_modes_prop) 2*num_modes]);
ticks= [1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop            1             floor(num_modes_prop/4)             floor(num_modes_prop/2)             floor(3*num_modes_prop/4)             num_modes_prop   num_modes_evanes];
yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes+1) (num_modes+floor(num_modes_prop/4)) (num_modes+floor(num_modes_prop/2)) (num_modes+floor(3*num_modes_prop/4)) (num_modes+num_modes_prop) 2*num_modes]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
title('$Phase(S)$','Interpreter','Latex')
set(gca,'FontSize',FontSizeVal)



% define a function to extract parts of Smatrix 
[mask_prop_prop,mask_evan_evan,mask_evan_prop,mask_prop_evan]=generate_selection_mask(init_data.num_modes,init_data.num_modes_prop);                       
            % mask_prop_prop : Incident mode : prop, scattered mode : evan
            % mask_evan_evan : Incident mode : prop, scattered mode : evan
            % mask_evan_prop : Incident mode : prop, scattered mode : evan
            % mask_prop_evan : Incident mode : evan, scattered mode : prop
% define a new partioned S matrix such that 
% Spartitioned= [Sprop,prop Sprop,evan;
%                Sevan,prop Sevan,evan] 
%------------------ Generalised unitarity relations -----------------------
Smatrix_prop_prop=reshape(Smatrix(mask_prop_prop==1),2*init_data.num_modes_prop,2*init_data.num_modes_prop);
Smatrix_evan_evan=reshape(Smatrix(mask_evan_evan==1),2*init_data.num_modes_evanes,2*init_data.num_modes_evanes);
Smatrix_evan_prop=reshape(Smatrix(mask_evan_prop==1),2*init_data.num_modes_evanes,2*init_data.num_modes_prop);
Smatrix_prop_evan=reshape(Smatrix(mask_prop_evan==1),2*init_data.num_modes_prop,2*init_data.num_modes_evanes);
[size(Smatrix_prop_prop) size(Smatrix_evan_evan) size(Smatrix_evan_prop) size(Smatrix_prop_evan)]
Spartitioned=[Smatrix_prop_prop Smatrix_prop_evan; Smatrix_evan_prop Smatrix_evan_evan];

figure('Position', [50 50 800 500],'color','W');
colormap jet
imagesc(abs(Spartitioned).^2)
axis equal tight
colorbar
set(gca,'FontSize',FontSizeVal) 
xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes_prop+floor(num_modes_prop/4)) (num_modes_prop+floor(num_modes_prop/2)) (num_modes_prop+floor(3*num_modes_prop/4)) (num_modes_prop+num_modes_prop) (2*num_modes_prop+num_modes_evanes) 2*num_modes]);
ticks= [1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop                floor(num_modes_prop/4)                  floor(num_modes_prop/2)                  floor(3*num_modes_prop/4)                  num_modes_prop                     num_modes_evanes  num_modes_evanes];
yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes_prop+floor(num_modes_prop/4)) (num_modes_prop+floor(num_modes_prop/2)) (num_modes_prop+floor(3*num_modes_prop/4)) (num_modes_prop+num_modes_prop) (2*num_modes_prop+num_modes_evanes) 2*num_modes]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
title('$|S_{parti}|^2$','Interpreter','Latex')
caxis([min(min((abs(Spartitioned).^2))) max(max((abs(Spartitioned).^2)))])
set(gca,'ColorScale','log')


Iprop=zeros(2*init_data.num_modes,2*init_data.num_modes);
Ievan=zeros(2*init_data.num_modes,2*init_data.num_modes);
Iprop(1:2*init_data.num_modes_prop,1:2*init_data.num_modes_prop)=eye(2*init_data.num_modes_prop,2*init_data.num_modes_prop);
Ievan(2*init_data.num_modes-2*init_data.num_modes_evanes+1:2*init_data.num_modes,2*init_data.num_modes-2*init_data.num_modes_evanes+1:2*init_data.num_modes)=eye(2*init_data.num_modes_evanes,2*init_data.num_modes_evanes);

figure('Position', [50 50 800 500],'color','W');
subplot(1,2,1)
spy(Iprop,'k')
xlabel('')
xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes_prop+floor(num_modes_prop/4)) (num_modes_prop+floor(num_modes_prop/2)) (num_modes_prop+floor(3*num_modes_prop/4)) (num_modes_prop+num_modes_prop) (2*num_modes_prop+num_modes_evanes) 2*num_modes]);
ticks= [1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop                floor(num_modes_prop/4)                  floor(num_modes_prop/2)                  floor(3*num_modes_prop/4)                  num_modes_prop                     num_modes_evanes  num_modes_evanes];
yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes_prop+floor(num_modes_prop/4)) (num_modes_prop+floor(num_modes_prop/2)) (num_modes_prop+floor(3*num_modes_prop/4)) (num_modes_prop+num_modes_prop) (2*num_modes_prop+num_modes_evanes) 2*num_modes]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
title('$I_{pr}$','Interpreter','Latex')
set(gca,'FontSize',FontSizeVal)

subplot(1,2,2)
spy(Ievan,'k')
xlabel('')
xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes_prop+floor(num_modes_prop/4)) (num_modes_prop+floor(num_modes_prop/2)) (num_modes_prop+floor(3*num_modes_prop/4)) (num_modes_prop+num_modes_prop) (2*num_modes_prop+num_modes_evanes) 2*num_modes]);
ticks= [1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop                floor(num_modes_prop/4)                  floor(num_modes_prop/2)                  floor(3*num_modes_prop/4)                  num_modes_prop                     num_modes_evanes  num_modes_evanes];
yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes_prop+floor(num_modes_prop/4)) (num_modes_prop+floor(num_modes_prop/2)) (num_modes_prop+floor(3*num_modes_prop/4)) (num_modes_prop+num_modes_prop) (2*num_modes_prop+num_modes_evanes) 2*num_modes]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
title('$I_{ev}$','Interpreter','Latex')
set(gca,'FontSize',FontSizeVal)


figure('Position', [50 50 800 500],'color','W');
imagesc(abs(Spartitioned'*Iprop*Spartitioned-1i.*(Spartitioned'*Ievan-Ievan*Spartitioned)));% Refer the paper or the thesis for the references.
title({'$Generalized~unitarity~relation$','$|S_{parti}^\dagger I_{pr} S_{parti} -i(S_{parti}^\dagger I_{ev}-I_{ev}S_{parti})|=I_{pr}$'},'Interpreter','Latex')
xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes_prop+floor(num_modes_prop/4)) (num_modes_prop+floor(num_modes_prop/2)) (num_modes_prop+floor(3*num_modes_prop/4)) (num_modes_prop+num_modes_prop) (2*num_modes_prop+num_modes_evanes) 2*num_modes]);
ticks= [1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop                floor(num_modes_prop/4)                  floor(num_modes_prop/2)                  floor(3*num_modes_prop/4)                  num_modes_prop                     num_modes_evanes  num_modes_evanes];
yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes_prop+floor(num_modes_prop/4)) (num_modes_prop+floor(num_modes_prop/2)) (num_modes_prop+floor(3*num_modes_prop/4)) (num_modes_prop+num_modes_prop) (2*num_modes_prop+num_modes_evanes) 2*num_modes]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
colormap(flipud(gray))
set(gca,'FontSize',FontSizeVal) 
colorbar

figure('Position', [50 50 800 500],'color','W');
imagesc(abs(Spartitioned'*Iprop*Spartitioned-1i.*(Spartitioned'*Ievan-Ievan*Spartitioned)-Iprop));% Refer the paper or the thesis for the references.
title({'$Generalized~unitarity~relation$','$|S_{parti}^\dagger I_{pr} S_{parti} -i(S_{parti}^\dagger I_{ev}-I_{ev}S_{parti})-I_{pr}|=0$'},'Interpreter','Latex')
xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes_prop+floor(num_modes_prop/4)) (num_modes_prop+floor(num_modes_prop/2)) (num_modes_prop+floor(3*num_modes_prop/4)) (num_modes_prop+num_modes_prop) (2*num_modes_prop+num_modes_evanes) 2*num_modes]);
ticks= [1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop                floor(num_modes_prop/4)                  floor(num_modes_prop/2)                  floor(3*num_modes_prop/4)                  num_modes_prop                     num_modes_evanes  num_modes_evanes];
yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes_prop+floor(num_modes_prop/4)) (num_modes_prop+floor(num_modes_prop/2)) (num_modes_prop+floor(3*num_modes_prop/4)) (num_modes_prop+num_modes_prop) (2*num_modes_prop+num_modes_evanes) 2*num_modes]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
set(gca,'FontSize',FontSizeVal) 
colorbar

%---------- Generalised reciprocity relation ------------------------------
figure('Position', [50 50 800 500],'color','W');
imagesc(abs(Smatrix_prop_prop-Smatrix_prop_prop.'));
% Refer the paper or the thesis for the references.
title({'$Generalized~reciprocity~relation$','$|S^{pr,pr}_{parti}-{(S^{pr,pr}_{parti})}^{T}|=0$'},'Interpreter','Latex')
xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes_prop+floor(num_modes_prop/4)) (num_modes_prop+floor(num_modes_prop/2)) (num_modes_prop+floor(3*num_modes_prop/4)) (num_modes_prop+num_modes_prop)]);
ticks= [1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop                floor(num_modes_prop/4)                  floor(num_modes_prop/2)                  floor(3*num_modes_prop/4)                  num_modes_prop ];
yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes_prop+floor(num_modes_prop/4)) (num_modes_prop+floor(num_modes_prop/2)) (num_modes_prop+floor(3*num_modes_prop/4)) (num_modes_prop+num_modes_prop)]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
axis equal tight
set(gca,'FontSize',FontSizeVal) 
colorbar

figure('Position', [50 50 800 500],'color','W');
imagesc(abs(Smatrix_prop_evan-1i.*Smatrix_evan_prop.'));% Refer the paper or the thesis for the references.
title({'$Generalized~reciprocity~relation$','$|S^{pr,ev}_{parti}-i{(S^{ev,pr}_{parti})}^{T}|=0$'},'Interpreter','Latex')
    yticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes_prop+floor(num_modes_prop/4)) (num_modes_prop+floor(num_modes_prop/2)) (num_modes_prop+floor(3*num_modes_prop/4)) (num_modes_prop+num_modes_prop)]);
ytickslab= [1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop                floor(num_modes_prop/4)                  floor(num_modes_prop/2)                  floor(3*num_modes_prop/4)                   num_modes_prop ];
    xticks([1 num_modes_evanes 2*num_modes_evanes]);
xtickslab=([1 num_modes_evanes   num_modes_evanes]);
xticklabels(xtickslab)
yticklabels(ytickslab)
xtickangle(90)
axis equal tight
set(gca,'FontSize',FontSizeVal) 
colorbar

figure('Position', [50 50 800 500],'color','W');
imagesc(abs(Smatrix_evan_prop+1i.*Smatrix_prop_evan.'));% Refer the paper or the thesis for the references.
title({'$Generalized~reciprocity~relation$','$|S^{ev,pr}_{parti}+i{(S^{pr,ev}_{parti})}^{T}|=0$'},'Interpreter','Latex')
    xticks([1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop (num_modes_prop+floor(num_modes_prop/4)) (num_modes_prop+floor(num_modes_prop/2)) (num_modes_prop+floor(3*num_modes_prop/4)) (num_modes_prop+num_modes_prop)]);
xtickslab= [1 floor(num_modes_prop/4) floor(num_modes_prop/2) floor(3*num_modes_prop/4) num_modes_prop                floor(num_modes_prop/4)                  floor(num_modes_prop/2)                  floor(3*num_modes_prop/4)                   num_modes_prop ];
    yticks([1 num_modes_evanes 2*num_modes_evanes]);
ytickslab=([1 num_modes_evanes   num_modes_evanes]);
yticklabels(ytickslab)
xticklabels(xtickslab)
ytickangle(90)
xtickangle(90)
axis equal tight
set(gca,'FontSize',FontSizeVal) 
colorbar

figure('Position', [50 50 800 500],'color','W');
imagesc(abs(Smatrix_evan_evan-Smatrix_evan_evan.'));
title({'$Generalized~reciprocity~relation$','$|S^{ev,ev}_{parti}-{(S^{ev,ev}_{parti})}^{T}|=0$'},'Interpreter','Latex')
xticks([1 num_modes_evanes  2*num_modes_evanes]);
ticks= [1 num_modes_evanes  num_modes_evanes];
yticks([1 num_modes_evanes  2*num_modes_evanes]);
xticklabels(ticks)
yticklabels(ticks)
xtickangle(90)
axis equal tight
set(gca,'FontSize',FontSizeVal) 
colorbar

%--- numerical example 1 involving evanescent & propagating modes ---------
mind=init_data.num_modes_prop-10;  % where m <= Mprop
nind=init_data.num_modes_prop+floor(init_data.num_modes_evanes/2); % n > Nprop                     
[S11(mind,nind) 1i.*S11(nind,mind) S11(mind,nind)-1i.*S11(nind,mind)]
[S22(mind,nind) 1i.*S22(nind,mind) S22(mind,nind)-1i.*S22(nind,mind)]
[S12(mind,nind) 1i.*S21(nind,mind) S12(mind,nind)-1i.*S21(nind,mind)]
[S21(mind,nind) 1i.*S12(nind,mind) S21(mind,nind)-1i.*S12(nind,mind)]

%----- numerical example 2 involving evanescent & propagating modes -------
nind=init_data.num_modes_prop-10;  % n <= Nprop
mind=init_data.num_modes_prop+floor(init_data.num_modes_evanes/2); % m> Mprop                     
[S11(mind,nind) -1i.*S11(nind,mind) S11(mind,nind)+1i.*S11(nind,mind)]
[S22(mind,nind) -1i.*S22(nind,mind) S22(mind,nind)+1i.*S22(nind,mind)]
[S12(mind,nind) -1i.*S21(nind,mind) S12(mind,nind)+1i.*S21(nind,mind)]
[S21(mind,nind) -1i.*S12(nind,mind) S21(mind,nind)+1i.*S12(nind,mind)]

end

