%Terzis Dimitrios AEM6101

%Post processor beam geranos 

%import initial data

%Check if nodes file exist
if exist('./AEM6101_Nodes_ger_beam.txt','file')==2
     nodes_beam=importdata('AEM6101_Nodes_ger_beam.txt',',');
else
     waitfor(msgbox('Nodes file not found!'))
     error('Nodes file not found')
end

%check if nodes new cord file exist
if exist('./AEM6101_NodesNewCord_ger_beam.txt','file')==2
     nodes_beam_newcord=importdata('AEM6101_NodesNewCord_ger_beam.txt',',');
else
     waitfor(msgbox('NodesNewCord file not found!'))
     error('NodesNewCord file not found')
end

%check if boundary conditions file exist
if exist('./AEM6101_BC_ger_beam.txt','file')==2
     bc_nodes_beam=importdata('AEM6101_BC_ger_beam.txt',',');
else
     waitfor(msgbox('BC file not found!'))
     error('BC file not found')
end


%check if forces file exist
if exist('./AEM6101_Forces_ger_beam.txt','file')==2
     force_beam=importdata('AEM6101_Forces_ger_beam.txt',',');
else
     waitfor(msgbox('Force file not found!'))
     error('Force file not found')
end

%check if elements file exist
if exist('./AEM6101_Elements_ger_beam.txt','file')==2
     elements_beam=importdata('AEM6101_Elements_ger_beam.txt',',');
else
     waitfor(msgbox('Elements file not found!'))
     error('Elements file not found')
end


%check if elements length file exist
if exist('./AEM6101_Elementslength_ger_beam.txt','file')==2
     elements_beam_length=importdata('AEM6101_Elementslength_ger_beam.txt',',');
else
     waitfor(msgbox('Elements length file not found!'))
     error('Elements length file not found')
end


%check if elements new length file exist
if exist('./AEM6101_ElementsNewlength_ger_beam.txt','file')==2
     elements_beam_newlength=importdata('AEM6101_ElementsNewlength_ger_beam.txt',',');
else
     waitfor(msgbox('Elements newlength file not found!'))
     error('Elements newlength file not found')
end

%check if Strains file exist
if exist('./AEM6101_Strains_ger_beam.txt','file')==2
     strain=importdata('AEM6101_Strains_ger_beam.txt',',');
else
     waitfor(msgbox('Strains file not found!'))
     error('Strains file not found')
end


%check if Displacement file exist
if exist('./AEM6101_Displacement_ger_beam.txt','file')==2
     U=importdata('AEM6101_Displacement_ger_beam.txt',',');
else
     waitfor(msgbox('Displacement file not found!'))
     error('Displacement file not found')
end

%check if Reactions file exist
if exist('./AEM6101_Reactions_ger_beam.txt','file')==2
     Reaction=importdata('AEM6101_Reactions_ger_beam.txt',',');
else
     waitfor(msgbox('Reactions file not found!'))
     error('Reactions file not found')
end


%check if Stresses file exist
if exist('./AEM6101_Stresses_ger_beam.txt','file')==2
     stress=importdata('AEM6101_Stresses_ger_beam.txt',',');
else
     waitfor(msgbox('Stresses file not found!'))
     error('Stresses file not found')
end

%cross sections
Ay=max(elements_beam(:,5));
Ad=min(elements_beam(:,5));

a=0;
s=10; %default scale factor
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
        if isnumeric(s)==1 %if input is numeric 
            a=1;
            s=answer;
        else
            waitfor(msgbox('Wrong Input! Use Number for Scale Factor')); %if wrong input while is continued
        end
    end
end

sc_Node_coords=nodes_beam(1:32,2:7)+s*U;

%-----------------PLOTS------------------
figure 
%----------first: 
%visualization of undeformed cube
subplot(2,1,1)
plot3(0,0,0,'o')
title('Undeformed Structure') %title of plots
xlabel('x-axis [mm]')   %specify units of measurement
ylabel('y-axis [mm]')   
zlabel('z-axis [mm]')
hold on
%create lines for each element
% line with 3 width for horizontal and vertical elements
% line with 1 width for diagonal elements

for i=1:size(elements_beam,1)
    if elements_beam(i,5)==Ay
        line([nodes_beam(elements_beam(i,2),2) nodes_beam(elements_beam(i,3),2)],[nodes_beam(elements_beam(i,2),3) nodes_beam(elements_beam(i,3),3)],[nodes_beam(elements_beam(i,2),4) nodes_beam(elements_beam(i,3),4)],'Color','black','LineWidth' ,3)
    else
        line([nodes_beam(elements_beam(i,2),2) nodes_beam(elements_beam(i,3),2)],[nodes_beam(elements_beam(i,2),3) nodes_beam(elements_beam(i,3),3)],[nodes_beam(elements_beam(i,2),4) nodes_beam(elements_beam(i,3),4)],'Color','black','LineWidth' ,1)
    end
end
%plot with circles each node
for i=1:size(nodes_beam,1)
    scatter3(nodes_beam(i,2),nodes_beam(i,3),nodes_beam(i,4),'filled','k')
end

view([-13 11])
camproj('orthographic')
axis equal
axis auto
hold off



%----------Second:
%visualization of Deformed Structure and Underformed together
subplot(2,1,2)
plot3(0,0,0,'o')
title('Undeformed - Deformed Structure') %title of plots
xlabel('x-axis [mm]')   %specify units of measurement
ylabel('y-axis [mm]')   
zlabel('z-axis [mm]')
hold on


%the undeformed
for i=1:size(elements_beam,1)
    if elements_beam(i,5)==Ay
        line([nodes_beam(elements_beam(i,2),2) nodes_beam(elements_beam(i,3),2)],[nodes_beam(elements_beam(i,2),3) nodes_beam(elements_beam(i,3),3)],[nodes_beam(elements_beam(i,2),4) nodes_beam(elements_beam(i,3),4)],'Color','black','LineStyle','--','LineWidth',1)
    else
        line([nodes_beam(elements_beam(i,2),2) nodes_beam(elements_beam(i,3),2)],[nodes_beam(elements_beam(i,2),3) nodes_beam(elements_beam(i,3),3)],[nodes_beam(elements_beam(i,2),4) nodes_beam(elements_beam(i,3),4)],'Color','black','LineStyle','--','LineWidth',1)
    end
end

%the deformed
for i=1:size(elements_beam,1)
    if elements_beam(i,5)==Ay
        line([sc_Node_coords(elements_beam(i,2),1) sc_Node_coords(elements_beam(i,3),1)],[sc_Node_coords(elements_beam(i,2),2) sc_Node_coords(elements_beam(i,3),2)],[sc_Node_coords(elements_beam(i,2),3) sc_Node_coords(elements_beam(i,3),3)],'Color','r','LineWidth' ,3)
    else
        line([sc_Node_coords(elements_beam(i,2),1) sc_Node_coords(elements_beam(i,3),1)],[sc_Node_coords(elements_beam(i,2),2) sc_Node_coords(elements_beam(i,3),2)],[sc_Node_coords(elements_beam(i,2),3) sc_Node_coords(elements_beam(i,3),3)],'Color','r','LineWidth' ,1)
    end
end

%plot each node
for i=1:(size(nodes_beam,1)-1)
    scatter3(sc_Node_coords(i,1),sc_Node_coords(i,2),sc_Node_coords(i,3),'filled','r')
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

%the undeformed
for i=1:size(elements_beam,1)
    if elements_beam(i,5)==Ay
        line([nodes_beam(elements_beam(i,2),2) nodes_beam(elements_beam(i,3),2)],[nodes_beam(elements_beam(i,2),3) nodes_beam(elements_beam(i,3),3)],[nodes_beam(elements_beam(i,2),4) nodes_beam(elements_beam(i,3),4)],'Color','black','LineStyle','--','LineWidth',3)
    else
        line([nodes_beam(elements_beam(i,2),2) nodes_beam(elements_beam(i,3),2)],[nodes_beam(elements_beam(i,2),3) nodes_beam(elements_beam(i,3),3)],[nodes_beam(elements_beam(i,2),4) nodes_beam(elements_beam(i,3),4)],'Color','black','LineStyle','--','LineWidth',1)
    end
end

%the deformed
ma = max(strain); 
mi = min(strain);
cmap = colormap(jet(256));
xx = linspace(0,1,256);
for i=1:size(elements_beam,1)
    se=(strain(i)-mi)/(ma-mi);
    if elements_beam(i,5)==Ay
        line([sc_Node_coords(elements_beam(i,2),1) sc_Node_coords(elements_beam(i,3),1)],[sc_Node_coords(elements_beam(i,2),2) sc_Node_coords(elements_beam(i,3),2)],[sc_Node_coords(elements_beam(i,2),3) sc_Node_coords(elements_beam(i,3),3)],'Color', interp1(xx, cmap, se),'LineWidth' ,3)
        %patch([sc_Node_coords(Elements.Nodes(i,1),1) sc_Node_coords(Elements.Nodes(i,2),1)],[sc_Node_coords(Elements.Nodes(i,1),2) sc_Node_coords(Elements.Nodes(i,2),2)],[sc_Node_coords(Elements.Nodes(i,1),3) sc_Node_coords(Elements.Nodes(i,2),3)],'EdgeColor','interp','LineWidth' ,3)
    else
        line([sc_Node_coords(elements_beam(i,2),1) sc_Node_coords(elements_beam(i,3),1)],[sc_Node_coords(elements_beam(i,2),2) sc_Node_coords(elements_beam(i,3),2)],[sc_Node_coords(elements_beam(i,2),3) sc_Node_coords(elements_beam(i,3),3)],'Color', interp1(xx, cmap, se),'LineWidth' ,2)
        %patch([sc_Node_coords(Elements.Nodes(i,1),1) sc_Node_coords(Elements.Nodes(i,2),1)],[sc_Node_coords(Elements.Nodes(i,1),2) sc_Node_coords(Elements.Nodes(i,2),2)],[sc_Node_coords(Elements.Nodes(i,1),3) sc_Node_coords(Elements.Nodes(i,2),3)],'EdgeColor','interp','LineWidth' ,1)
    end
end

%plot each node
for i=1:(size(nodes_beam,1)-1)
    scatter3(sc_Node_coords(i,1),sc_Node_coords(i,2),sc_Node_coords(i,3),'filled','b')
end


caxis([mi ma]); h = colorbar;
ylabel(h, 'Axial Strain')
camproj('orthographic')
view([-13 11])
axis auto
axis equal
hold off


%--------------Fourth:
%the deformed structure with Axis Stresses
subplot(2,1,2)
plot3(0,0,0,'o')
title('Deformed Structure with Axial Stresses') %title of plots
xlabel('x-axis [mm]')   %specify units of measurement
ylabel('y-axis [mm]')   
zlabel('z-axis [mm]')
hold on

%the undeformed
for i=1:size(elements_beam,1)
    if elements_beam(i,5)==Ay
        line([nodes_beam(elements_beam(i,2),2) nodes_beam(elements_beam(i,3),2)],[nodes_beam(elements_beam(i,2),3) nodes_beam(elements_beam(i,3),3)],[nodes_beam(elements_beam(i,2),4) nodes_beam(elements_beam(i,3),4)],'Color','black','LineStyle','--','LineWidth',3)
    else
        line([nodes_beam(elements_beam(i,2),2) nodes_beam(elements_beam(i,3),2)],[nodes_beam(elements_beam(i,2),3) nodes_beam(elements_beam(i,3),3)],[nodes_beam(elements_beam(i,2),4) nodes_beam(elements_beam(i,3),4)],'Color','black','LineStyle','--','LineWidth',1)
    end
end

%the deformed
ma = max(stress); 
mi = min(stress);
cmap = colormap(jet(256));
xx = linspace(0,1,256);
for i=1:size(elements_beam,1)
    se=(stress(i)-mi)/(ma-mi);
    if elements_beam(i,5)==Ay
        line([sc_Node_coords(elements_beam(i,2),1) sc_Node_coords(elements_beam(i,3),1)],[sc_Node_coords(elements_beam(i,2),2) sc_Node_coords(elements_beam(i,3),2)],[sc_Node_coords(elements_beam(i,2),3) sc_Node_coords(elements_beam(i,3),3)],'Color', interp1(xx, cmap, se),'LineWidth' ,3)
        %patch([sc_Node_coords(Elements.Nodes(i,1),1) sc_Node_coords(Elements.Nodes(i,2),1)],[sc_Node_coords(Elements.Nodes(i,1),2) sc_Node_coords(Elements.Nodes(i,2),2)],[sc_Node_coords(Elements.Nodes(i,1),3) sc_Node_coords(Elements.Nodes(i,2),3)],'EdgeColor','interp','LineWidth' ,3)
    else
        line([sc_Node_coords(elements_beam(i,2),1) sc_Node_coords(elements_beam(i,3),1)],[sc_Node_coords(elements_beam(i,2),2) sc_Node_coords(elements_beam(i,3),2)],[sc_Node_coords(elements_beam(i,2),3) sc_Node_coords(elements_beam(i,3),3)],'Color', interp1(xx, cmap, se),'LineWidth' ,2)
        %patch([sc_Node_coords(Elements.Nodes(i,1),1) sc_Node_coords(Elements.Nodes(i,2),1)],[sc_Node_coords(Elements.Nodes(i,1),2) sc_Node_coords(Elements.Nodes(i,2),2)],[sc_Node_coords(Elements.Nodes(i,1),3) sc_Node_coords(Elements.Nodes(i,2),3)],'EdgeColor','interp','LineWidth' ,1)
    end
end

%plot each node
for i=1:(size(nodes_beam,1)-1)
    scatter3(sc_Node_coords(i,1),sc_Node_coords(i,2),sc_Node_coords(i,3),'filled','b')
end


caxis([mi ma]);
h = colorbar;
ylabel(h, 'Axial Stress')
camproj('orthographic')
view([-13 11])
axis auto
axis equal
hold off



figure
%-----------------Fifth
%the deformed structure with displacements x
subplot(2,2,1)
plot3(0,0,0,'o')
title('Deformed Structure on Ux') %title of plots
xlabel('x-axis [mm]')   %specify units of measurement
ylabel('y-axis [mm]')   
zlabel('z-axis [mm]')
hold on


%the undeformed
for i=1:size(elements_beam,1)
    if elements_beam(i,5)==Ay
        line([nodes_beam(elements_beam(i,2),2) nodes_beam(elements_beam(i,3),2)],[nodes_beam(elements_beam(i,2),3) nodes_beam(elements_beam(i,3),3)],[nodes_beam(elements_beam(i,2),4) nodes_beam(elements_beam(i,3),4)],'Color','black','LineStyle','--','LineWidth',1)
    else
        line([nodes_beam(elements_beam(i,2),2) nodes_beam(elements_beam(i,3),2)],[nodes_beam(elements_beam(i,2),3) nodes_beam(elements_beam(i,3),3)],[nodes_beam(elements_beam(i,2),4) nodes_beam(elements_beam(i,3),4)],'Color','black','LineStyle','--','LineWidth',1)
    end
end


%the deformed
Ux=U(:,1);
ma=max(Ux);
mi=min(Ux);
cmap = colormap(jet(1000));
xx = linspace(0,1,1000);
for i=1:size(elements_beam,1)
    if Ux(elements_beam(i,2))>Ux(elements_beam(i,3))
        cmin=(Ux(elements_beam(i,3))-mi)/(ma-mi); %Nodes coefficient with lowest displacement 
        cmax=(Ux(elements_beam(i,2))-mi)/(ma -mi); %Nodes coefficient with higher displacement
        kx=linspace(sc_Node_coords(elements_beam(i,3),1),sc_Node_coords(elements_beam(i,2),1),200); %from lowest to highest
        ky=linspace(sc_Node_coords(elements_beam(i,3),2),sc_Node_coords(elements_beam(i,2),2),200);
        kz=linspace(sc_Node_coords(elements_beam(i,3),3),sc_Node_coords(elements_beam(i,2),3),200);
        a=linspace(cmin,cmax,200);
        for j=2:200
            if elements_beam(i,5)==Ay
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,3) 
            else
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,1)
            end
        end
    else
        cmin=(Ux(elements_beam(i,2))-mi)/(ma-mi); %Nodes coefficient with lowest displacement 
        cmax=(Ux(elements_beam(i,3))-mi)/(ma -mi); %Nodes coefficient with higher displacement
        kx=linspace(sc_Node_coords(elements_beam(i,2),1),sc_Node_coords(elements_beam(i,3),1),200); %from lowest to highest
        ky=linspace(sc_Node_coords(elements_beam(i,2),2),sc_Node_coords(elements_beam(i,3),2),200);
        kz=linspace(sc_Node_coords(elements_beam(i,2),3),sc_Node_coords(elements_beam(i,3),3),200);
        a=linspace(cmin,cmax,200);
        for j=2:200
            if elements_beam(i,5)==Ay
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,3) 
            else
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,1)
            end
        end
    end
end


caxis([mi ma]);
h = colorbar;
ylabel(h, 'Ux Displacements')
camproj('orthographic')
view([-13 11])
axis auto
axis equal
hold off


%-----------------Sixth
%the deformed structure with displacements y
subplot(2,2,2)
plot3(0,0,0,'o')
title('Deformed Structure on Uy') %title of plots
xlabel('x-axis [mm]')   %specify units of measurement
ylabel('y-axis [mm]')   
zlabel('z-axis [mm]')
hold on

%the undeformed
for i=1:size(elements_beam,1)
    if elements_beam(i,5)==Ay
        line([nodes_beam(elements_beam(i,2),2) nodes_beam(elements_beam(i,3),2)],[nodes_beam(elements_beam(i,2),3) nodes_beam(elements_beam(i,3),3)],[nodes_beam(elements_beam(i,2),4) nodes_beam(elements_beam(i,3),4)],'Color','black','LineStyle','--','LineWidth',1)
    else
        line([nodes_beam(elements_beam(i,2),2) nodes_beam(elements_beam(i,3),2)],[nodes_beam(elements_beam(i,2),3) nodes_beam(elements_beam(i,3),3)],[nodes_beam(elements_beam(i,2),4) nodes_beam(elements_beam(i,3),4)],'Color','black','LineStyle','--','LineWidth',1)
    end
end

%the deformed
Uy=U(:,2);
ma=max(Uy);
mi=min(Uy);
cmap = colormap(jet(1000));
xx = linspace(0,1,1000);
for i=1:size(elements_beam,1)
    if Uy(elements_beam(i,2))>Uy(elements_beam(i,3))
        cmin=(Uy(elements_beam(i,3))-mi)/(ma-mi); %Nodes coefficient with lowest displacement 
        cmax=(Uy(elements_beam(i,2))-mi)/(ma -mi); %Nodes coefficient with higher displacement
        kx=linspace(sc_Node_coords(elements_beam(i,3),1),sc_Node_coords(elements_beam(i,2),1),200); %from lowest to highest
        ky=linspace(sc_Node_coords(elements_beam(i,3),2),sc_Node_coords(elements_beam(i,2),2),200);
        kz=linspace(sc_Node_coords(elements_beam(i,3),3),sc_Node_coords(elements_beam(i,2),3),200);
        a=linspace(cmin,cmax,200);
        for j=2:200
            if elements_beam(i,5)==Ay
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,3) 
            else
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,1)
            end
        end
    else
        cmin=(Uy(elements_beam(i,2))-mi)/(ma-mi); %Nodes coefficient with lowest displacement 
        cmax=(Uy(elements_beam(i,3))-mi)/(ma -mi); %Nodes coefficient with higher displacement
        kx=linspace(sc_Node_coords(elements_beam(i,2),1),sc_Node_coords(elements_beam(i,3),1),200); %from lowest to highest
        ky=linspace(sc_Node_coords(elements_beam(i,2),2),sc_Node_coords(elements_beam(i,3),2),200);
        kz=linspace(sc_Node_coords(elements_beam(i,2),3),sc_Node_coords(elements_beam(i,3),3),200);
        a=linspace(cmin,cmax,200);
        for j=2:200
            if elements_beam(i,5)==Ay
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,3) 
            else
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,1)
            end
        end
    end
end



caxis([mi ma]);
h = colorbar;
ylabel(h, 'Uy Displacements')
camproj('orthographic')
view([-13 11])
axis auto
axis equal
hold off


%---------------Seventh
%the deformed structure with displacements z
subplot(2,2,3)
plot3(0,0,0,'o')
title('Deformed Structure on Uz') %title of plots
xlabel('x-axis [mm]')   %specify units of measurement
ylabel('y-axis [mm]')   
zlabel('z-axis [mm]')
hold on

%the undeformed
for i=1:size(elements_beam,1)
    if elements_beam(i,5)==Ay
        line([nodes_beam(elements_beam(i,2),2) nodes_beam(elements_beam(i,3),2)],[nodes_beam(elements_beam(i,2),3) nodes_beam(elements_beam(i,3),3)],[nodes_beam(elements_beam(i,2),4) nodes_beam(elements_beam(i,3),4)],'Color','black','LineStyle','--','LineWidth',1)
    else
        line([nodes_beam(elements_beam(i,2),2) nodes_beam(elements_beam(i,3),2)],[nodes_beam(elements_beam(i,2),3) nodes_beam(elements_beam(i,3),3)],[nodes_beam(elements_beam(i,2),4) nodes_beam(elements_beam(i,3),4)],'Color','black','LineStyle','--','LineWidth',1)
    end
end

%the deformed
Uz=U(:,3);
ma=max(Uz);
mi=min(Uz);
cmap = colormap(jet(1000));
xx = linspace(0,1,1000);
for i=1:size(elements_beam,1)
    if Uz(elements_beam(i,2))>Uz(elements_beam(i,3))
        cmin=(Uz(elements_beam(i,3))-mi)/(ma-mi); %Nodes coefficient with lowest displacement 
        cmax=(Uz(elements_beam(i,2))-mi)/(ma -mi); %Nodes coefficient with higher displacement
        kx=linspace(sc_Node_coords(elements_beam(i,3),1),sc_Node_coords(elements_beam(i,2),1),200); %from lowest to highest
        ky=linspace(sc_Node_coords(elements_beam(i,3),2),sc_Node_coords(elements_beam(i,2),2),200);
        kz=linspace(sc_Node_coords(elements_beam(i,3),3),sc_Node_coords(elements_beam(i,2),3),200);
        a=linspace(cmin,cmax,200);
        for j=2:200
            if elements_beam(i,5)==Ay
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,3) 
            else
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,1)
            end
        end
    else
        cmin=(Uz(elements_beam(i,2))-mi)/(ma-mi); %Nodes coefficient with lowest displacement 
        cmax=(Uz(elements_beam(i,3))-mi)/(ma -mi); %Nodes coefficient with higher displacement
        kx=linspace(sc_Node_coords(elements_beam(i,2),1),sc_Node_coords(elements_beam(i,3),1),200); %from lowest to highest
        ky=linspace(sc_Node_coords(elements_beam(i,2),2),sc_Node_coords(elements_beam(i,3),2),200);
        kz=linspace(sc_Node_coords(elements_beam(i,2),3),sc_Node_coords(elements_beam(i,3),3),200);
        a=linspace(cmin,cmax,200);
        for j=2:200
            if elements_beam(i,5)==Ay
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,3) 
            else
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,1)
            end
        end
    end
end


caxis([mi ma]);
h = colorbar;
ylabel(h, 'Uz Displacements')
camproj('orthographic')
view([-13 11])
axis auto
axis equal
hold off

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
for i=1:size(elements_beam,1)
    if elements_beam(i,5)==Ay
        line([nodes_beam(elements_beam(i,2),2) nodes_beam(elements_beam(i,3),2)],[nodes_beam(elements_beam(i,2),3) nodes_beam(elements_beam(i,3),3)],[nodes_beam(elements_beam(i,2),4) nodes_beam(elements_beam(i,3),4)],'Color','black','LineStyle','--','LineWidth',1)
    else
        line([nodes_beam(elements_beam(i,2),2) nodes_beam(elements_beam(i,3),2)],[nodes_beam(elements_beam(i,2),3) nodes_beam(elements_beam(i,3),3)],[nodes_beam(elements_beam(i,2),4) nodes_beam(elements_beam(i,3),4)],'Color','black','LineStyle','--','LineWidth',1)
    end
end

%the deformed
Ut=sqrt(Ux.^2+Uy.^2+Uz.^2);
ma=max(Ut);
mi=min(Ut);
cmap = colormap(jet(1000));
xx = linspace(0,1,1000);

for i=1:size(elements_beam,1)
    if Ut(elements_beam(i,2))>Ut(elements_beam(i,3))
        cmin=(Ut(elements_beam(i,3))-mi)/(ma-mi); %Nodes coefficient with lowest displacement 
        cmax=(Ut(elements_beam(i,2))-mi)/(ma -mi); %Nodes coefficient with higher displacement
        kx=linspace(sc_Node_coords(elements_beam(i,3),1),sc_Node_coords(elements_beam(i,2),1),200); %from lowest to highest
        ky=linspace(sc_Node_coords(elements_beam(i,3),2),sc_Node_coords(elements_beam(i,2),2),200);
        kz=linspace(sc_Node_coords(elements_beam(i,3),3),sc_Node_coords(elements_beam(i,2),3),200);
        a=linspace(cmin,cmax,200);
        for j=2:200
            if elements_beam(i,5)==Ay
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,3) 
            else
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,1)
            end
        end
    else
        cmin=(Ut(elements_beam(i,2))-mi)/(ma-mi); %Nodes coefficient with lowest displacement 
        cmax=(Ut(elements_beam(i,3))-mi)/(ma -mi); %Nodes coefficient with higher displacement
        kx=linspace(sc_Node_coords(elements_beam(i,2),1),sc_Node_coords(elements_beam(i,3),1),200); %from lowest to highest
        ky=linspace(sc_Node_coords(elements_beam(i,2),2),sc_Node_coords(elements_beam(i,3),2),200);
        kz=linspace(sc_Node_coords(elements_beam(i,2),3),sc_Node_coords(elements_beam(i,3),3),200);
        a=linspace(cmin,cmax,200);
        for j=2:200
            if elements_beam(i,5)==Ay
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,3) 
            else
                line([kx(j-1) kx(j)],[ky(j-1) ky(j)],[kz(j-1) kz(j)],'Color', interp1(xx, cmap, a(j)),'LineWidth' ,1)
            end
        end
    end
end



caxis([mi ma]);
h = colorbar;
ylabel(h, 'Utotal Displacements')
camproj('orthographic')
view([-13 11])
axis auto
axis equal
hold off

