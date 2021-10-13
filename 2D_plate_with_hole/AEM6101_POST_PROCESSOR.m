%% Post Processor
%Terzis Dimitrios 
%Aem 6101

%% Data  Input

clear all

%%
%-----------------------------Input_TXT_Data-------------------------------
%Check if nodes file exist
if exist('./AEM6101_Nodes_Plate.txt','file')==2
     node=importdata('AEM6101_Nodes_Plate.txt',',',1);
else
     waitfor(msgbox('Nodes file not found!'))
     error('Nodes file not found')
end

%check if elements file exist
if exist('./AEM6101_Elements_Plate.txt','file')==2
     element=importdata('AEM6101_Elements_Plate.txt',',',1);
else
     waitfor(msgbox('Elements file not found!'))
     error('Elements file not found')
end

%check if data file exist
if exist('./AEM6101_data_Plate.txt','file')==2
     data=importdata('AEM6101_data_Plate.txt',',',1);
else
     waitfor(msgbox('data file not found!'))
     error('data file not found')
end

%check if Strains file exist
if exist('./AEM6101_Strains_Plate.txt','file')==2
     strain=importdata('AEM6101_Strains_Plate.txt',',',1);
else
     waitfor(msgbox('Strains file not found!'))
     error('Strains file not found')
end

%check if Displacements file exist
if exist('./AEM6101_Displacement_Plate.txt','file')==2
     disp=importdata('AEM6101_Displacement_Plate.txt',',',1);
else
     waitfor(msgbox('Displacement file not found!'))
     error('Displacement file not found')
end


%check if Stresses file exist
if exist('./AEM6101_Stresses_Plate.txt','file') == 2
     stress = importdata('AEM6101_Stresses_Plate.txt',',',1);
else
     waitfor(msgbox('Stresses file not found!'))
     error('Stresses file not found')
end


%%
%initialize nodes struct and import data node to nodes struct
Nodes = struct('id', [], 'coords', [], 'newcoords', [], 'BCs', [], 'F', [], 'inner', [], 'perimeter', [],'sx', [],'sy', [],'txy',[],'s1',[],'s2',[],'phi',[],'vM',[]) ;
Nodes.id = node.data(:, 1); %nodes id
Nodes.coords = node.data(:, 2 : 3); %nodes coordinates
Nodes.newcoords = node.data(:, 8 : 9); %nodes new coordinates
Nodes.BCs = node.data(:, 4 : 5); % nodes boundary conditions
Nodes.F = node.data(:, 6 : 7); %nodes force
Nodes.Flags = node.data(:, 10); %nodes flags

%initialize elements struct and import data element to elements struct
Elements = struct('id',[],'Nodes',[],'A',[],'CG',[],'t',[]);
Elements.id = element.data(:,1); %elements id
Elements.Nodes = element.data(:,2:4); %nodes of its element
Elements.A = element.data(:,5); %elements cross-section
Elements.CG = element.data(:,6:7); %elements center of mass
Elements.t = element.data(:,8); %elements width
Elements.CGdef = element.data(:,9:10); %new Center of mass
Elements
E = data.data(1); %elasticity
v = data.data(2); %Poisson
G = data.data(3); %Torsion Stiffnes
a = data.data(4); %half length
b = data.data(5); %half width
x_nodes = data.data(6); %x nodes
y_nodes = data.data(7); %y nodes
d = data.data(10);   %diameter
t = data.data(11);   %thickness
stress_type = data.data(8); %1 for tensile, 2 for bending
Load_nom = data.data(9);    %Nomina load value      


%initialize stresses struct and import data element to elements struct
Strain = struct('ex',[],'ey',[],'gxy',[]);
Strain.ex = strain.data(:,1);
Strain.ey = strain.data(:,2);
Strain.gxy = strain.data(:,3);

%initialize stresses struct and import data element to elements struct
Stress = struct('sx',[],'sy',[],'txy',[],'s1',[],'s2',[],'phi',[],'vM',[]);
Stress.sx = stress.data(:,1);
Stress.sy = stress.data(:,2);
Stress.txy = stress.data(:,3);
Stress.s1 = stress.data(:,4);
Stress.s2 = stress.data(:,5);
Stress.phi = stress.data(:,6);
Stress.vM = stress.data(:,7);

%initialiaze struct for displacements

%input needed
x = 0;
%while x == 0
    %scale = str2double(inputdlg('Plese define scale factor'));
    %if isnumeric(scale) == 1
        %x = 1;
    %else
        %waitfor(msgbox('Make sure your answer is numerric'));
    %end
%end
scale = 500;
U = struct('x',[],'y',[],'total',[]);
U.x = disp.data(:, 1) ;
U.y = disp.data(:, 2) ;
for i = 1 : size(U.x, 1)
    U.total(i, 1) = norm([U.x(i) U.y(i)]) ;
end


clear data element node stress strain disp 


%% Extrapolation to the nodes
%first we need to seperate the perimeter nodes from the rest
temp = Nodes.Flags;
deik = find(Nodes.Flags == 2);
%Nodes.Flags([1, y_nodes, y_nodes*x_nodes+ 1, y_nodes * x_nodes + 1 - y_nodes, deik(1), deik(ceil(x_nodes/2 + 1))  ]) = 6;
Nodes.Flags(find(Nodes.Flags == 2)) = 0;
Nodes.Flags(find(Nodes.Flags == 4)) = 0;



Nodes.inner = struct('id', [], 'x', [], 'y', [], 'msx', [], 'msy', [], 'mtxy', [], 'ms1', [], 'ms2', [], 'mvM', []); 
Nodes.inner.id = (Nodes.id(Nodes.Flags == 0));
Nodes.inner.x = Nodes.coords(Nodes.Flags == 0, 1);
Nodes.inner.y = Nodes.coords(Nodes.Flags == 0, 2);



for i = 1:size(Nodes.inner.id)
    
    [Row, col] = find(Elements.Nodes(:, 1) == Nodes.inner.id(i)); %For Node #, find elements in which its present.
    Nodes.inner.msx(i, 1) = mean(Stress.sx(Row)); % Sx of that Node, is the mean of its surrounding elements.
    Nodes.inner.msy(i, 1) = mean(Stress.sy(Row));
    Nodes.inner.mtxy(i, 1) = mean(Stress.txy(Row));
    Nodes.inner.ms1(i, 1) = mean(Stress.s1(Row));
    Nodes.inner.ms2(i, 1) = mean(Stress.s2(Row));
    Nodes.inner.mvM(i, 1) = mean(Stress.vM(Row));
end




Nodes.perimeter = struct('id', [], 'x', [], 'y', [], 'msx', [], 'msy', [], 'mtxy', [], 'ms1', [], 'ms2', [], 'mvM', []);
Nodes.perimeter.id = (Nodes.id(Nodes.Flags ~= 0));
Nodes.perimeter.x = Nodes.coords(Nodes.Flags ~= 0, 1);
Nodes.perimeter.y = Nodes.coords(Nodes.Flags ~= 0, 2);



for i = 1:size(Nodes.perimeter.id, 1)
    [Row, col] = find(Elements.Nodes(:, 1) == Nodes.perimeter.id(i)); %For Node #, find elements in which its present.
    if numel(find(Elements.Nodes(:, 1) == Nodes.perimeter.id(i))) == 0
        [Row, col] = find(Elements.Nodes(:, 3) == Nodes.perimeter.id(i));
    end
    Nodes.perimeter.msx(i, 1)  = max(abs(Stress.sx(Row, 1)) .*sum(Stress.sx(Row, 1)) / abs(sum(Stress.sx(Row, 1)))); % Sx of that Node, is the mean of its surrounding elements.
    % to *sum(stres(row ktl ktl einai wste na exoyme to swsto prosimo sthn
    % tash mas. alliws to max pairnei thn pio megalh arithmitika kai alliw
    %s ola ginontai thetika. to sum einai giati sta kentrika elements,
    %kapoia exoyn thetiko, kapoia arnitiko prosimo tashs. Ara vriskoume to
    %prosimo apo to athrisma gia kalyteri proseggisi.
    Nodes.perimeter.msy(i, 1) = max(Stress.sy(Row, 1).*sum(Stress.sy(Row, 1)) / abs(sum(Stress.sy(Row, 1))));
    Nodes.perimeter.mtxy(i, 1) = max(Stress.txy(Row, 1).*sum(Stress.txy(Row, 1)) / abs(sum(Stress.txy(Row, 1))));
    Nodes.perimeter.ms1(i, 1) = max(Stress.s1(Row, 1).*sum(Stress.s1(Row, 1)) / abs(sum(Stress.s1(Row, 1))));
    Nodes.perimeter.ms2(i, 1) = max(Stress.s2(Row, 1).*sum(Stress.s2(Row, 1)) / abs(sum(Stress.s2(Row, 1))));
    Nodes.perimeter.mvM(i, 1) = max(Stress.vM(Row, 1)); % aytos einai mono thetikos 
end
Nodes.Flags = temp;


%% taseis Sta NODES
Nodes.sx = zeros(size(Nodes.id,1),1);
Nodes.sx(Nodes.perimeter.id(:,1)) = Nodes.perimeter.msx(:, 1);
Nodes.sx(Nodes.inner.id(:,1)) = Nodes.inner.msx(:, 1);

Nodes.sy = zeros(size(Nodes.id,1),1);
Nodes.sy(Nodes.perimeter.id(:,1)) = Nodes.perimeter.msy(:, 1);
Nodes.sy(Nodes.inner.id(:,1)) = Nodes.inner.msy(:, 1);

Nodes.txy = zeros(size(Nodes.id,1),1);
Nodes.txy(Nodes.perimeter.id(:,1)) = Nodes.perimeter.mtxy(:, 1);
Nodes.txy(Nodes.inner.id(:,1)) = Nodes.inner.mtxy(:, 1);

Nodes.s1 = zeros(size(Nodes.id,1),1);
Nodes.s1(Nodes.perimeter.id(:,1)) = Nodes.perimeter.ms1(:, 1);
Nodes.s1(Nodes.inner.id(:,1)) = Nodes.inner.ms1(:, 1);

Nodes.s2 = zeros(size(Nodes.id,1),1);
Nodes.s2(Nodes.perimeter.id(:,1)) = Nodes.perimeter.ms2(:, 1);
Nodes.s2(Nodes.inner.id(:,1))= Nodes.inner.ms2(:, 1);

Nodes.vM = zeros(size(Nodes.id,1),1);
Nodes.vM(Nodes.perimeter.id(:,1)) = Nodes.perimeter.mvM(:, 1);
Nodes.vM(Nodes.inner.id(:,1)) = Nodes.inner.mvM(:, 1);

Nodes.phi = atan(2.*Nodes.txy./(Nodes.sx-Nodes.sy))./2;

Nodes.sx(find(isnan(Nodes.sx)))=0;
%%
total_n = size(Nodes.id);
total_el = size(Elements.id);
v = 0.33;


sc_Nodes_coords(:, :) = [Nodes.coords(:, 1) + scale.*U.x(:) Nodes.coords(:, 2)+scale.*U.y(:)]; 
for i = 1 : size(Elements.Nodes, 1)
    
    sc_Elements_CGs(i, 1) = sum(sc_Nodes_coords(Elements.Nodes(i, 1:3), 1)) / 3;
    sc_Elements_CGs(i, 2) = sum(sc_Nodes_coords(Elements.Nodes(i, 1:3), 2)) / 3;
end


%% Notch Factors
if stress_type == 1
    Kt = 3 - 3.14*(d / 100) + 3.667*(d / 100)^ 2 - 1.527 * (d / 100)^2
    Snom = Load_nom / ((100-d)*4)
    Smax_FEM = Nodes.sx(1)
    Notch_FEM = Smax_FEM / Snom
    error_perc = abs(Notch_FEM - Kt) / Kt * 100
    waitfor(msgbox(sprintf('Calculated Notch factor: %d', Notch_FEM)))
    waitfor(msgbox(sprintf('Error percent: %d', error_perc)))
    
   
else
    Kt = 2;
    Snom_nohole = Load_nom * 12 * 50 / (100^3 * 4);
    Snom = d * (6 * Load_nom) / ((100^3 - d^3)*4)
    Smax_FEM = Nodes.sx(1)
    Notch_hole = Smax_FEM / Snom
    error_perc = abs((Notch_hole - Kt) / Kt) * 100
    waitfor(msgbox(sprintf('Calculated Notch factor: %d', Notch_hole)))
    waitfor(msgbox(sprintf('Error percent: %d', error_perc)))
    
    Kt_edge = 2*d/50;
    Snom_edge = (6*Load_nom*50) / ((100^3 - d^3)*4)
    Smax_FEM_edge = max(Nodes.sx(find(Nodes.Flags==5)))
    Notch_edge = Smax_FEM_edge / Snom_edge

end


%% Plots
y = 0;
while y == 0
    [indx, y] = listdlg('PromptString', 'Select plot target', 'SelectionMode', 'single', 'ListString',{'Plot on nodes', 'Plot on elements'});
end
x = 0; 
while x == 0
    if indx == 1 %plot on nodes
        [request, x] = listdlg('PromptString','Select Result to plot','SelectionMode','single','ListString',{'Displacement on x', 'Displacement on y', 'Total Displacement', 'sx', 'sy', 'txy', 's1', 's2', 'vonMises'});
    else        %plot on elements
        [request, x] = listdlg('PromptString','Select Result to plot','SelectionMode','single','ListString',{'sx', 'sy', 'txy', 'ex', 'ey', 'gxy', 's1', 's2', 'vonMises'});
    end
end

%plot on elements needs colours and calculated forces on elements
%plot on nodes needs nodes_colours and calculated forces on nodes

colours = [Stress.sx Stress.sy Stress.txy Strain.ex Strain.ey Strain.gxy Stress.s1 Stress.s2 Stress.vM];
titles = ["sx", "sy", "txy", "ex", "ey" , "gxy", "s1", "s2", "vonMises"];

nodes_titles = ["Displacement on x", "Displacement on y", "Total Displacement", "sx", "sy", "txy", "s1", "s2", "vonMises"];
nodes_colours = [U.x U.y U.total Nodes.sx Nodes.sy Nodes.txy Nodes.s1 Nodes.s2 Nodes.vM];

x_labels = "mm";
y_labels = "mm";



if indx == 2
    Plot_maker(request, titles, Nodes.coords, sc_Nodes_coords, Elements.Nodes, x_labels, y_labels, colours(:, request), indx, Elements.CG, sc_Elements_CGs, Stress.phi)
else
    Plot_maker(request, nodes_titles, Nodes.coords, sc_Nodes_coords, Elements.Nodes, x_labels, y_labels, nodes_colours(:, request), indx, Elements.CGdef, sc_Elements_CGs, Stress.phi)
end


%%

function[] = Plot_maker(request, titles, Nodes, NewNodes, Elements, x_labels, y_labels, colours, indx, old_cgs, new_cgs, Elements_phi, Nodes_phi)

if indx == 2 
    figure
    patch('Faces',Elements,'Vertices',Nodes,'FaceVertexCData',[colours],'FaceColor','none')
    axis equal
    title('Undeformed');
    xlabel = x_labels;
    ylabel = y_labels;
    if request == 7 || request == 8
        hold on
        Quiver_Maker(Elements_phi, colours, old_cgs )
    end


    figure
    patch('Faces', Elements, 'Vertices', NewNodes, 'FaceVertexCData', [colours], 'FaceColor', 'flat')
    axis equal
    dx = colorbar;
    dx.Ticks = [linspace(min(colours), max(colours), 13)];
    caxis([min(colours) max(colours)])
    axis equal
    title(titles(request))
    %title = (answer);
    xlabel = x_labels;
    ylabel = y_labels;
    if request == 7 || request == 8
        hold on
        Quiver_Maker(Elements_phi, colours, new_cgs )
    end
    
    

else
    figure
    patch('Faces',Elements,'Vertices',Nodes,'FaceVertexCData',[colours],'FaceColor','none')
    axis equal
    title('Undeformed');
    xlabel = x_labels;
    ylabel = y_labels;
    if request == 8 || request == 9
        hold on
        Quiver_Maker(Nodes_phi, colours, old_cgs)
    end


    figure
    patch('Faces', Elements, 'Vertices', NewNodes, 'FaceVertexCData', [colours], 'FaceColor', 'interp')
    axis equal
    dx = colorbar;
    dx.Ticks = [linspace(min(colours), max(colours), 13)];
    caxis([min(colours) max(colours)])
    axis equal
    title(titles(request))
    %title = (answer);
    xlabel = x_labels;
    ylabel = y_labels;
    if request == 8 || request == 9
        hold on
        Quiver_Maker(Nodes_phi, colours, new_cgs)
    end
end

    
end


%%
function[] = Quiver_Maker(phi, variable, CGs)
hold on
q = quiver(CGs(:, 1), CGs(:, 2), cos(phi) .* variable(:), sin(phi) .* variable(:));
c = q.Color;
q.Color = 'black';
end

