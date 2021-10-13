%geranos post processor
%Terzis Dimitrios AEM6101

%Check if nodes file exist
if exist('./AEM6101_Nodes_ger.txt','file')==2
     nodes=importdata('AEM6101_Nodes_ger.txt',',');
else
     waitfor(msgbox('Nodes file not found!'))
     error('Nodes file not found')
end

%check if nodes new cord file exist
if exist('./AEM6101_NodesNewCord_ger.txt','file')==2
     nodes_newcord=importdata('AEM6101_NodesNewCord_ger.txt',',');
else
     waitfor(msgbox('NodesNewCord file not found!'))
     error('NodesNewCord file not found')
end

%check if boundary conditions file exist
if exist('./AEM6101_BC_ger.txt','file')==2
     bc_nodes=importdata('AEM6101_BC_ger.txt',',');
else
     waitfor(msgbox('BC file not found!'))
     error('BC file not found')
end

%check if forces file exist
if exist('./AEM6101_Forces_ger.txt','file')==2
     force=importdata('AEM6101_Forces_ger.txt',',');
else
     waitfor(msgbox('Force file not found!'))
     error('Force file not found')
end

%check if elements file exist
if exist('./AEM6101_Elements_ger.txt','file')==2
     elements=importdata('AEM6101_Elements_ger.txt',',');
else
     waitfor(msgbox('Elements file not found!'))
     error('Elements file not found')
end

%check if elements length file exist
if exist('./AEM6101_Elementslength_ger.txt','file')==2
     elements_length=importdata('AEM6101_Elementslength_ger.txt',',');
else
     waitfor(msgbox('Elements length file not found!'))
     error('Elements length file not found')
end

%check if elements new length file exist
if exist('./AEM6101_ElementsNewlength_ger.txt','file')==2
     elements_newlength=importdata('AEM6101_Elementslength_ger.txt',',');
else
     waitfor(msgbox('Elements newlength file not found!'))
     error('Elements newlength file not found')
end


%check if Strains file exist
if exist('./AEM6101_Strains_ger.txt','file')==2
     strain=importdata('AEM6101_Strains_ger.txt',',');
else
     waitfor(msgbox('Strains file not found!'))
     error('Strains file not found')
end


%check if Displacement file exist
if exist('./AEM6101_Displacement_ger.txt','file')==2
     U=importdata('AEM6101_Displacement_ger.txt',',');
else
     waitfor(msgbox('Displacement file not found!'))
     error('Displacement file not found')
end

%check if Reactions file exist
if exist('./AEM6101_Reactions_ger.txt','file')==2
     Reaction=importdata('AEM6101_Reactions_ger.txt',',');
else
     waitfor(msgbox('Reactions file not found!'))
     error('Reactions file not found')
end

%check if Stresses file exist
if exist('./AEM6101_Stresses_ger.txt','file')==2
     stress=importdata('AEM6101_Stresses_ger.txt',',');
else
     waitfor(msgbox('Stresses file not found!'))
     error('Stresses file not found')
end


nodes=nodes.data;
nodes_newcord=nodes_newcord.data;
bc_nodes=bc_nodes.data;
force=force.data;
elements=elements.data;
elements_length=elements_length.data;
elements_newlength=elements_newlength.data;
strain=strain.data;
U=U.data;
Reaction=Reaction.data;
stress=stress.data;
a=0;
scale=10; %default scale factor
while a~=1
    [indx,tf]=listdlg('ListString',{'Input Scale Factor','Use default scale factor (10)'},'SelectionMode','single');
    if tf==1
        a=1;
    else
        waitfor(msgbox('A selection is needed!'));
    end
end
if indx==1
    a=0; % variable for while
    dims=[1 35];
    while a~=1
        answer=inputdlg({'Enter scale factor'},'Input',dims,{'10'});
        answer=str2double(answer);
        if isnumeric(scale)==1 %if input is numeric 
            a=1;
            scale=answer;
        else
            waitfor(msgbox('Wrong Input! Use Number for Scale Factor')); %if wrong input while is continued
        end
    end
end

Ay=max(elements(:,4));
scaled_Node_coords=nodes(:,2:4)+scale*U;

%-----------------PLOTS------------------
figure
%----------first: 
%visualization of undeformed cube
subplot(1,2,1)
plot3(0,0,0,'o')
title('Undeformed Structure') %title of plots
xlabel('x-axis [mm]')   %specify units of measurement
ylabel('y-axis [mm]')   
zlabel('z-axis [mm]')
hold on
%create lines for each element
% line with 3 width for horizontal and vertical elements


for i=1:size(elements,1)
        line([nodes(elements(i,2),2) nodes(elements(i,3),2)],[nodes(elements(i,2),3) nodes(elements(i,3),3)],[nodes(elements(i,2),4) nodes(elements(i,3),4)],'Color','black','LineWidth' ,3)
end
%plot with circles each node
for i=1:size(nodes,1)
    scatter3(nodes(i,2),nodes(i,3),nodes(i,4),'filled','k')
end
view([-13 11])
camproj('orthographic')
axis equal
axis auto
hold off



%----------Second:
%visualization of Deformed Structure and Underformed together
subplot(1,2,2)
plot3(0,0,0,'o')
title('Undeformed - Deformed Structure') %title of plots
xlabel('x-axis [mm]')   %specify units of measurement
ylabel('y-axis [mm]')   
zlabel('z-axis [mm]')
hold on


%undeformed
for i=1:size(elements,1)
    line([nodes(elements(i,2),2) nodes(elements(i,3),2)],[nodes(elements(i,2),3) nodes(elements(i,3),3)],[nodes(elements(i,2),4) nodes(elements(i,3),4)],'Color','black','LineStyle','--','LineWidth',3)
end

%deformed
for i=1:size(elements,1)
    line(scaled_Node_coords(elements(i,2:3),1),scaled_Node_coords(elements(i,2:3),2),scaled_Node_coords(elements(i,2:3),3),'Color','r','LineWidth' ,1)
end

%plot each node
for i=1:size(nodes,1)
    scatter3(scaled_Node_coords(i,1),scaled_Node_coords(i,2),scaled_Node_coords(i,3),'filled','r')
end

view([-13 11])
camproj('orthographic')
axis equal
axis auto
hold off




figure
%------------Third: 
%the deformed structure with Axis Strain
subplot(2,1,1)
plot3(0,0,0,'o')
title('Deformed Structure with Axial Strains') %title of plots
xlabel('x-axis [mm]')   %specify units of measurement
ylabel('y-axis [mm]')   
zlabel('z-axis [mm]')
hold on

%undeformed
for i=1:size(elements,1)
        line(nodes(elements(i,2:3),2),nodes(elements(i,2:3),3),nodes(elements(i,2:3),4),'Color','black','LineStyle','--','LineWidth',1)
end

%deformed
max_strain = max(strain); 
min_strain = min(strain);
cmap = colormap(jet(1500));
xx = linspace(0,1,1500);
for i=1:size(elements,1)
    strain_local=(strain(i)-min_strain)/(max_strain-min_strain);
    line(scaled_Node_coords(elements(i,2:3),1),scaled_Node_coords(elements(i,2:3),2),scaled_Node_coords(elements(i,2:3),3),'Color', interp1(xx, cmap, strain_local),'LineWidth' ,3)      
end

%plot each node
for i=1:size(nodes,1)
    scatter3(scaled_Node_coords(i,1),scaled_Node_coords(i,2),scaled_Node_coords(i,3),'filled','r')
end
caxis([min_strain max_strain]); h = colorbar;
ylabel(h, 'Axial Strain')
axis auto
axis equal
hold off


%--------------Fourth:
%deformed structure with Axis Stresses
subplot(2,1,2)
plot3(0,0,0,'o')
title('Deformed Structure with Axial Stresses') %title of plots
xlabel('x-axis [mm]')   %specify units of measurement
ylabel('y-axis [mm]')   
zlabel('z-axis [mm]')
hold on

%undeformed
for i=1:size(elements,1)
        line(nodes(elements(i,2:3),2),nodes(elements(i,2:3),3),nodes(elements(i,2:3),4),'Color','black','LineStyle','--','LineWidth',1)
end

%deformed
ma = max(stress); 
mi = min(stress);
cmap = colormap(jet(1500));
xx = linspace(0,1,1500);
for i=1:size(elements,1)
    se=(stress(i)-mi)/(ma-mi);
    line([scaled_Node_coords(elements(i,2),1) scaled_Node_coords(elements(i,3),1)],[scaled_Node_coords(elements(i,2),2) scaled_Node_coords(elements(i,3),2)],[scaled_Node_coords(elements(i,2),3) scaled_Node_coords(elements(i,3),3)],'Color', interp1(xx, cmap, se),'LineWidth' ,3)
end

%plot each node
for i=1:size(nodes,1)
    scatter3(scaled_Node_coords(i,1),scaled_Node_coords(i,2),scaled_Node_coords(i,3),'filled','r')
end

caxis([mi ma]);
h = colorbar;
axis auto
axis equal
hold off


figure
%-----------------Fifth
%deformed structure with displacements x
Ux=U(:,1);
if sum(Ux(:))==0
else
    subplot(1,3,1)
    plot3(0,0,0,'o')
    title('Deformed Structure on Ux') %title of plots
    xlabel('x-axis [mm]')   %specify units of measurement
    ylabel('y-axis [mm]')   
    zlabel('z-axis [mm]')
    hold on

%undeformed
    for i=1:size(elements,1)
        line(nodes(elements(i,2:3),2),nodes(elements(i,2:3),3),nodes(elements(i,2:3),4),'Color','black','LineStyle','--','LineWidth',1)
    end

%deformed
    
    max_deffx=max(Ux);
    min_deffx=min(Ux);
    cmap = colormap(jet(1000));
    xx = linspace(0,1,1000);
    for i=1:size(elements,1)
        cmin=(Ux(elements(i,3))-min_deffx)/(max_deffx-min_deffx); %Nodes coefficient with lowest displacement 
        cmax=(Ux(elements(i,2))-min_deffx)/(max_deffx-min_deffx); %Nodes coefficient with higher displacement
        a=linspace(cmin,cmax,200);
        if Ux(elements(i,2))>Ux(elements(i,3))
            kx=linspace(scaled_Node_coords(elements(i,3),1),scaled_Node_coords(elements(i,2),1),200); %from lowest to highest
            ky=linspace(scaled_Node_coords(elements(i,3),2),scaled_Node_coords(elements(i,2),2),200);
            kz=linspace(scaled_Node_coords(elements(i,3),3),scaled_Node_coords(elements(i,2),3),200);
            for j=2:200
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,3)
            end
        else
            kx=linspace(scaled_Node_coords(elements(i,2),1),scaled_Node_coords(elements(i,3),1),200); %from lowest to highest
            ky=linspace(scaled_Node_coords(elements(i,2),2),scaled_Node_coords(elements(i,3),2),200);
            kz=linspace(scaled_Node_coords(elements(i,2),3),scaled_Node_coords(elements(i,3),3),200);
            for j=2:200
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,3) 
            end
        end
    end


    caxis([min_deffx max_deffx]);
    h = colorbar;
    ylabel(h, 'Ux Displacements')
    axis auto
    axis equal
    hold off
end

%-----------------Sixth
%the deformed structure with displacements y

Uy=U(:,2);
if sum(Uy(:))
    subplot(1,3,2)
    plot3(0,0,0,'o')
    title('Deformed Structure on Uy') %title of plots
    xlabel('x-axis [mm]')   %specify units of measurement
    ylabel('y-axis [mm]')   
    zlabel('z-axis [mm]')
    hold on

%undeformed
    for i=1:size(elements,1)
        line(nodes(elements(i,2:3),2),nodes(elements(i,2:3),3),nodes(elements(i,2:3),4),'Color','black','LineStyle','--','LineWidth',1)
    end

%deformed
    
    max_deffy=max(Uy);
    min_deffy=min(Uy);
    cmap = colormap(jet(1000));
    xx = linspace(0,1,1000);
    for i=1:size(elements,1)
        cmin=(Uy(elements(i,3))-min_deffy)/(max_deffy-min_deffy); %Nodes coefficient with lowest displacement 
        cmax=(Uy(elements(i,2))-min_deffy)/(max_deffy-min_deffy); %Nodes coefficient with higher displacement
        a=linspace(cmin,cmax,200);
        if Uy(elements(i,2))>Uy(elements(i,3))
            kx=linspace(scaled_Node_coords(elements(i,3),1),scaled_Node_coords(elements(i,2),1),200); %from lowest to highest
            ky=linspace(scaled_Node_coords(elements(i,3),2),scaled_Node_coords(elements(i,2),2),200);
            kz=linspace(scaled_Node_coords(elements(i,3),3),scaled_Node_coords(elements(i,2),3),200);
            for j=2:200
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,3)
            end
        else
            kx=linspace(scaled_Node_coords(elements(i,3),1),scaled_Node_coords(elements(i,2),1),200); %from lowest to highest
            ky=linspace(scaled_Node_coords(elements(i,3),2),scaled_Node_coords(elements(i,2),2),200);
            kz=linspace(scaled_Node_coords(elements(i,3),3),scaled_Node_coords(elements(i,2),3),200);
            for j=2:200
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,3) 
            
            end
        end
    end


    caxis([min_deffy max_deffy]);
    h = colorbar;
	ylabel(h, 'Uy Displacements')
    axis auto
    axis equal
    hold off
end







%-----------------Seventh
%the deformed structure with displacements z

Uz=U(:,3);
if sum(Uz(:))==0
else
    subplot(1,3,3)
    plot3(0,0,0,'o')
    title('Deformed Structure on Uz') %title of plots
    xlabel('x-axis [mm]')   %specify units of measurement
    ylabel('y-axis [mm]')   
    zlabel('z-axis [mm]')
    hold on

%undeformed
    for i=1:size(elements,1)
        line(nodes(elements(i,2:3),2),nodes(elements(i,2:3),3),nodes(elements(i,2:3),4),'Color','black','LineStyle','--','LineWidth',1)
    end


%deformed
    max_deffz=max(Uz);
    min_deffz=min(Uz);
    cmap = colormap(jet(1000));
    xx = linspace(0,1,1000);
    for i=1:size(elements,1)
        cmin=(Uz(elements(i,3))-min_deffz)/(max_deffz-min_deffz); %Nodes coefficient with lowest displacement 
        cmax=(Uz(elements(i,2))-min_deffz)/(max_deffz-min_deffz); %Nodes coefficient with higher displacement
        a=linspace(cmin,cmax,200);
        if Uz(elements(i,2))>Uz(elements(i,3))
            kx=linspace(scaled_Node_coords(elements(i,3),1),scaled_Node_coords(elements(i,2),1),200); %from lowest to highest
            ky=linspace(scaled_Node_coords(elements(i,3),2),scaled_Node_coords(elements(i,2),2),200);
            kz=linspace(scaled_Node_coords(elements(i,3),3),scaled_Node_coords(elements(i,2),3),200);
            for j=2:200
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,3)
            end
        else
            kx=linspace(scaled_Node_coords(elements(i,3),1),scaled_Node_coords(elements(i,2),1),200); %from lowest to highest
            ky=linspace(scaled_Node_coords(elements(i,3),2),scaled_Node_coords(elements(i,2),2),200);
            kz=linspace(scaled_Node_coords(elements(i,3),3),scaled_Node_coords(elements(i,2),3),200);
            for j=2:200
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,3) 
            
            end
        end
    end


    caxis([min_deffz max_deffz]);
    h = colorbar;
    ylabel(h, 'Uz Displacements')
    axis auto
    axis equal
    hold off
end




%---------------Eighth
%the deformed structure with total displacement
subplot(2,2,4)
plot3(0,0,0,'o')
title('Deformed Structure on Utotal') %title of plots
xlabel('x-axis [mm]')   %specify units of measurement
ylabel('y-axis [mm]')   
zlabel('z-axis [mm]')
hold on

%the undeformed
for i=1:size(elements_ger,1)
    if elements_ger(i,4)==Ay
        line([nodes_ger(elements_ger(i,2),2) nodes_ger(elements_ger(i,3),2)],[nodes_ger(elements_ger(i,2),3) nodes_ger(elements_ger(i,3),3)],[nodes_ger(elements_ger(i,2),4) nodes_ger(elements_ger(i,3),4)],'Color','black','LineStyle','--','LineWidth',1)
    else
        line([nodes_ger(elements_ger(i,2),2) nodes_ger(elements_ger(i,3),2)],[nodes_ger(elements_ger(i,2),3) nodes_ger(elements_ger(i,3),3)],[nodes_ger(elements_ger(i,2),4) nodes_ger(elements_ger(i,3),4)],'Color','black','LineStyle','--','LineWidth',1)
    end
end

%the deformed
Ut=sqrt(Ux.^2+Uy.^2+Uz.^2);
ma=max(Ut);
mi=min(Ut);
cmap = colormap(jet(1000));
xx = linspace(0,1,1000);
for i=1:size(elements_ger,1)
    if Ut(elements_ger(i,2))>Ut(elements_ger(i,3))
        cmin=(Ut(elements_ger(i,3))-mi)/(ma-mi); %Nodes coefficient with lowest displacement 
        cmax=(Ut(elements_ger(i,2))-mi)/(ma -mi); %Nodes coefficient with higher displacement
        kx=linspace(scaled_Node_coords(elements_ger(i,3),1),scaled_Node_coords(elements_ger(i,2),1),200); %from lowest to highest
        ky=linspace(scaled_Node_coords(elements_ger(i,3),2),scaled_Node_coords(elements_ger(i,2),2),200);
        kz=linspace(scaled_Node_coords(elements_ger(i,3),3),scaled_Node_coords(elements_ger(i,2),3),200);
        a=linspace(cmin,cmax,200);
        for j=2:200
            if elements_ger(i,4)==Ay
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,3) 
            else
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,1)
            end
        end
    else
        cmin=(Ut(elements_ger(i,2))-mi)/(ma-mi); %Nodes coefficient with lowest displacement 
        cmax=(Ut(elements_ger(i,3))-mi)/(ma -mi); %Nodes coefficient with higher displacement
        kx=linspace(scaled_Node_coords(elements_ger(i,2),1),scaled_Node_coords(elements_ger(i,3),1),200); %from lowest to highest
        ky=linspace(scaled_Node_coords(elements_ger(i,2),2),scaled_Node_coords(elements_ger(i,3),2),200);
        kz=linspace(scaled_Node_coords(elements_ger(i,2),3),scaled_Node_coords(elements_ger(i,3),3),200);
        a=linspace(cmin,cmax,200);
        for j=2:200
            if elements_ger(i,4)==Ay
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,3) 
            else
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,1)
            end
        end
    end
     caxis([min_deffz max_deffz]);
    h = colorbar;
    ylabel(h, 'Uôtotal Displacements')
    axis auto
    axis equal
    hold off
end

