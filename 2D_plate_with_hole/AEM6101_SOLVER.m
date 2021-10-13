%Terzis Dimitrios
%AEM 6101
%solver
clear all
close all


%%
%-----------------------------Input_TXT_Data-------------------------------
%Check if nodes file exist
if exist('./AEM6101_Nodes_Plate.txt', 'file') == 2
     node = importdata('AEM6101_Nodes_Plate.txt', ',' ,1);
else
     waitfor(msgbox('Nodes file not found!'))
     error('Nodes file not found')
end

%check if elements file exist
if exist('./AEM6101_Elements_Plate.txt', 'file') == 2
     element = importdata('AEM6101_Elements_Plate.txt', ',',  1);
else
     waitfor(msgbox('Elements file not found!'))
     error('Elements file not found')
end

%check if data file exist
if exist('./AEM6101_data_Plate.txt', 'file') == 2
     data = importdata('AEM6101_data_Plate.txt', ',', 1);
else
     waitfor(msgbox('Data file not found!'))
     error('Data file not found')
end
%%
%initialize nodes struct and import data node to nodes struct
Nodes = struct('id',[],'coords',[],'newcoords',[],'BCs',[],'F',[]);
Nodes.id = node.data(:, 1); %nodes id
Nodes.coords = node.data(:, 2 : 3); %nodes coordinates
Nodes.BCs = node.data(:, 4 : 5); % nodes boundary conditions
Nodes.F = node.data(:, 6 : 7); %nodes force
Nodes.Flag = node.data(:, 8); %nodes flags

%initialize elements struct and import data element to elements struct
Elements = struct('id', [], 'Nodes', [],'A', [], 'CG', [], 't', [], 'CGdef', [], 'e', [], 's', [], 'ms', [], 'vM', []);
Elements.id = element.data(:, 1); %elements id
Elements.Nodes = element.data(:, 2 : 4); %nodes of its element
Elements.A = element.data(:, 5); %elements cross-section
Elements.CG = element.data(:, 6 : 7); %elements center of mass
Elements.t = element.data(:, 8); %elements width

E=data.data(1); %elasticity
v=data.data(2); %Poisson
G=data.data(3); %Torsion Stiffnes
a=data.data(4); %half length
b=data.data(5); %half width


clear data element node
%K_stiff = zeros(2*size(Nodes.id,1));

%%
%Stiffness


[B D K_stiff] = AEM6101_stiffness(Nodes.coords, Elements.Nodes, E, v, 4);




K = K_stiff;

for i = 1 : length(Nodes.id)
    if Nodes.BCs(i, 1) == 1
        K(2 * i - 1, 2 * i - 1) = K(2 * i - 1, 2 * i - 1) + 10^7;
    end
    if Nodes.BCs(i, 2) == 1
        K(2 * i, 2 * i) = K(2 * i,2 * i) + 10^7;
    end
end

F=reshape(Nodes.F', [], 1);
%% Solver 

%Displacements

u = K \ F;
U = reshape(u, 2, [])';
Nodes.newcoords = Nodes.coords + U;

%Reactions
Rs = K_stiff * u;
Reactions = reshape(Rs, 2, [])';



%%
for i = 1 : size(Elements.id, 1)
    b1 = Nodes.coords(Elements.Nodes(i, 2), 2) - Nodes.coords(Elements.Nodes(i, 3), 2);
    b2 = Nodes.coords(Elements.Nodes(i, 3), 2) - Nodes.coords(Elements.Nodes(i, 1), 2);
    b3 = Nodes.coords(Elements.Nodes(i, 1), 2) - Nodes.coords(Elements.Nodes(i, 2), 2);
    
    c1 = Nodes.coords(Elements.Nodes(i, 3), 1) - Nodes.coords(Elements.Nodes(i, 2), 1);
    c2 = Nodes.coords(Elements.Nodes(i, 1), 1) - Nodes.coords(Elements.Nodes(i, 3), 1);
    c3 = Nodes.coords(Elements.Nodes(i, 2), 1) - Nodes.coords(Elements.Nodes(i, 1), 1);
    
    B = 1 / (2 * Elements.A(i, 1)) * [b1 0 b2 0 b3 0; 0 c1 0 c2 0 c3; c1 b1 c2 b2 c3 b3];
    D = E / (1 - v^2) * [1 v 0 ; v 1 0; 0 0 (1 - v)/2];
    
    Elements.CGdef(i, 1) = sum(Nodes.newcoords(Elements.Nodes(i, :), 1)) / 3; 
    Elements.CGdef(i, 2) = sum(Nodes.newcoords(Elements.Nodes(i, :), 2)) / 3;
    
    Elements.e(i,:) = (B * [U(Elements.Nodes(i, 1), :) U(Elements.Nodes(i, 2), :) U(Elements.Nodes(i, 3), :)]')';
    Elements.s(i, :) = (D*(Elements.e(i,:))')';
    Elements.ms(i, 1) = (Elements.s(i, 1) + Elements.s(i,2)) / 2 + sqrt(((Elements.s(i, 1)-Elements.s(i, 2))/2)^2 + Elements.s(i, 3)^2);
    Elements.ms(i, 2) = (Elements.s(i, 1) + Elements.s(i, 2)) / 2 - sqrt(((Elements.s(i, 1)-Elements.s(i, 2))/2)^2 + Elements.s(i, 3)^2);
    Elements.ms(i, 3) = atan(2*Elements.s(i, 3) / (Elements.s(i, 1) - Elements.s(i, 2))) / 2;
    Elements.vM(i, 1) = sqrt(Elements.s(i,1)^2 - Elements.s(i, 1) * Elements.s(i, 2) + Elements.s(i, 2)^2 + 3*Elements.s(i, 3)^2);    
end


%%
fID=fopen('AEM6101_Nodes_Plate.txt','w');
fprintf(fID,'%3s,%13.5s,%13.5s,%5s,%5s,%13.8s,%13.8s,%13.5s,%13.5s,%3s \n','id','x','y','BCx','BCy','Fx','Fy','xnew','ynew', 'Flags');
for i=1:size(Nodes.id)
    fprintf(fID,'%3d,%13.5d,%13.5d,%5d,%5d,%13.8d,%13.8d,%13.5d,%13.5d,%3d \n',Nodes.id(i),Nodes.coords(i,1),Nodes.coords(i,2),Nodes.BCs(i,1),Nodes.BCs(i,2),Nodes.F(i,1),Nodes.F(i,2),Nodes.newcoords(i,1),Nodes.newcoords(i,2), Nodes.Flag(i));
end
fclose(fID);



%create txt file for elements
fID=fopen('AEM6101_Elements_Plate.txt','w');
fprintf(fID,'%5s,%4s,%4s,%4s,%8.2s,%8.4s,%8.4s,%4s,%8.4s,%8.4s \n','id','N1','N2','N3','A','CGx','CGy','t','CGdefx','CGdefy');
for i=1:size(Elements.id)
    fprintf(fID,'%5d,%4d,%4d,%4d,%8.4d,%8.4d,%8.4d,%4d,%8.4d,%8.4d \n',Elements.id(i,1), Elements.Nodes(i,1), Elements.Nodes(i,2), Elements.Nodes(i,3), Elements.A(i,1), Elements.CG(i,1), Elements.CG(i,2), Elements.t(i,1), Elements.CGdef(i,1), Elements.CGdef(i,2));
end
fclose(fID);


%create txt file for Displacements
fID = fopen('AEM6101_Displacement_Plate.txt', 'w');
fprintf(fID,'%11.4s,%11.4s \n', 'Ux', 'Uy');
for i = 1:size(U, 1)
    fprintf(fID, '%11.4d, %11.4d \n', U(i, 1), U(i, 2));
end

%create Reactions txt file
fID=fopen('AEM6101_Reactions_Plate.txt', 'w');
fprintf(fID,'%15.2s, %15.2s \n', 'Rx', 'Ry');
for i = 1 : size(Reactions,1)
    fprintf(fID,'%15.2d,%15.2d \n', Reactions(i, 1), Reactions(i, 2));
end

%create Strains txt file
fID = fopen('AEM6101_Strains_Plate.txt', 'w');
fprintf(fID, '%12.8s,%12.8s,%12.8s \n', 'ex', 'ey', 'gxy');
for i = 1 : size(Elements.e, 1)
    fprintf(fID, '%12.4d,%12.4d,%12.4d \n', Elements.e(i, 1), Elements.e(i, 2), Elements.e(i, 3));
end

%create txt file for Stresses
fID=fopen('AEM6101_Stresses_Plate.txt','w');
fprintf(fID,'%12.7s,%12.7s,%12.7s,%12.7s,%12.7s,%12.7s,%12.7s \n','sx','sy','txy','s1','s2','phi','VonMises');
for i=1:size(Elements.s,1)
    fprintf(fID,'%12.4d,%12.4d,%12.4d,%12.4d,%12.4d,%12.4d,%12.4d \n',Elements.s(i,1),Elements.s(i,2),Elements.s(i,3),Elements.ms(i,1),Elements.ms(i,2),Elements.ms(i,3),Elements.vM(i,1));
end
%%





function [B_global,D_global,STIFFNESS] = AEM6101_stiffness(NodesCoordinates,ElementNodes,E,v,t)

Nodof = 2;      % Degrees Of Freedom
NumberOfNodes = length(NodesCoordinates);
Stiffness=zeros(NumberOfNodes*Nodof);       % Initialize Stiffness matrice
B=[];                                       % Initialize B matrice
D=[];                                       % Initialize D matrice
for i=1:size(ElementNodes,1)
    node_1 = ElementNodes(i,1);
    node_2 = ElementNodes(i,2);
    node_3 = ElementNodes(i,3);
    x1 = NodesCoordinates(node_1,1); x2 = NodesCoordinates(node_2,1); x3 = NodesCoordinates(node_3,1);
    y1 = NodesCoordinates(node_1,2); y2 = NodesCoordinates(node_2,2); y3 = NodesCoordinates(node_3,2);
    A = (1/2)*det([1 x1 y1;1 x2 y2;1 x3 y3]);
    b1 = y2-y3;
    b2 = y3-y1;
    b3 = y1-y2;
    c1 = x3-x2;
    c2 = x1-x3;
    c3 = x2-x1;
    Be = (1/(2*A))*[b1 0 b2 0 b3 0;0 c1 0 c2 0 c3;c1 b1 c2 b2 c3 b3];            % Element B matrix
    B = [B;Be];                                                                  % Global B matrix
    De = (E/(1-(v^2)))*[1 v 0;v 1 0;0 0 (1-v)/2];                                % Element D matrix
    D = [D;De];                                                                  % Global D matrix
    index = [2*node_1-1 2*node_1 2*node_2-1 2*node_2 2*node_3-1 2*node_3];       % Getting the indice for the stiffness matrix
    Stiffness(index,index) = Stiffness(index,index)+ Be'*De*Be*t*A;              % Computing stiffness matrix
end
STIFFNESS = Stiffness;
B_global = B;
D_global = D;


end


