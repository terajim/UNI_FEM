clear all
close all
%%Terzis Dimitrios
%AEM 6101
%Pre - Processor
%% perioxi input 

%AEM input
q = 0;
while q == 0
    AEM = inputdlg('Write down your AEM', 'AEM', [1 15], "6101");
    AEM = str2double(AEM);
    if isnumeric(AEM) == 1
        q = 1;
    else
        waitfor(msgdlg('Wrong input, please try again'))
    end
end
digits=zeros(1 ,4);
digits(1 ,1) = floor(AEM / 1000);
digits(1 ,2) = floor((AEM - digits(1 ,1) * 1000) /100);
digits(1 ,3) = floor((AEM - digits(1 ,1) * 1000 - digits(1 ,2) * 100) / 10);
digits(1 ,4) = AEM - digits(1 ,1)*1000 - digits(1 ,2)*100 - digits(1 ,3) * 10;
d = 25 + 20 * (10 * digits(4) + digits(3)) / 100;
t = 4;

%Stoixeia Ylikou
E = 210000;
v = 0.33;
G = E / (2 * (1 + v));
d = 57;
lenx = 85; % half length x dir
leny = 50; % half length y dir



dims = [1 45];
q = 0;
while q == 0
    [list, x] = listdlg('PromptString','Please select a meshing mode','SelectionMode','single','ListString',{'Fast Meshing', 'Absolute Freedom'});
    if x == 1
        q = 1;
    else
        errordlg('You have to select a mode.')
        error('Pre-Processor terminated due to user input')
    end

end


if list == 1 
    q = 0;
    while q == 0
        answer = inputdlg({'x-axis nodes at the quarter plate', 'y-axis nodes at the quarter plate',...
            'x-axis spacing', 'y-axis spacing', 'hole spacing'},'Input region', dims, {'17', '10', '1', '1', '1'});
        answer = str2double(answer);
        if sum(isnan(answer)) == 0 && size(find(answer), 1) == size(answer, 1)
            q = 1;
        else
            waitfor(msgbox('Please fill in all blanks with non zero values'))
        end
    end
    x_nodes = answer(1);
    y_nodes = answer(2);
    x_spacing = answer(3);
    y_spacing = answer(4);
    hole_spacing = answer(5);
    
    flat_bottom_spacing = x_spacing;
    right_spacing = y_spacing;
    left_spacing = y_spacing;
    flat_top_spacing = x_spacing;
    percentage = 0.5;
else
    q = 0;
    while q == 0
        answer = inputdlg({'x-axis nodes at the quarter plate', 'y-axis nodes at the quarter plate',...
            'flat bottom region spacing', 'right region spacing','left region spacing', 'flat upper region spacing',...
            'hole spacing', '% of x-axis nodes at hole'},'Input region', dims, {'17', '10', '1', '1', '1', '1', '1', '0.5'});
        answer = str2double(answer);
        if sum(isnan(answer)) == 0 && size(find(answer), 1) == size(answer, 1)
            q = 1;
        else
            waitfor(msgbox('Please fill in all blanks with non zero values'))
        end
    end
    x_nodes = answer(1);
    y_nodes = answer(2);
    flat_bottom_spacing = answer(3);
    right_spacing = answer(4);
    left_spacing = answer(5);
    flat_top_spacing = answer(6);
    hole_spacing = answer(7);
    percentage = answer(8);
    
    y_spacing = (right_spacing + left_spacing) / 2;
    x_spacing = (flat_bottom_spacing + flat_top_spacing) / 2;    
end



%%
%katamoirasmos ARI8MOU nodes se mikos kai platos

n_do_flat = floor((1 -percentage) * x_nodes);       %katw flat meros

n_l_circle = floor(0.6 * (x_nodes - n_do_flat)) + 1;   %katw meros tou kuklou
n_u_circle = ceil(0.4 * (x_nodes - n_do_flat) ) + 1;    %panw meros kuklou

n_le_flat = y_nodes;       %left side nodes
n_r_flat = y_nodes;       %right side nodes 

n_up_flat = x_nodes;      %uper side nodes

nodes_number = [1 n_u_circle n_l_circle n_le_flat n_r_flat n_up_flat n_do_flat];




n_total_perimeter = (2 * (x_nodes + y_nodes));
arithmisi = transpose(1:n_total_perimeter);
nodes = zeros(n_total_perimeter, 3);
nodes(1:end, 1) = arithmisi(1:end, 1);
a8roistis = 0;
gwnies = [pi/2, pi/3, 0];


epilogi = listdlg('PromptString','Please select a meshing mode','SelectionMode','single','ListString',["Plate with hole", "Plate Without"]);
switch epilogi
    case 1
        perimeter = [d/2 0; 85 0; 85 50; 0 50; 0 d/2; d/2*cosd(60) d/2*sind(60)];     
    for i = 1 : 2
        a8roistis = a8roistis + nodes_number(i+1)-1;
        auto = dspace(gwnies(i), gwnies(i + 1), nodes_number(i + 1), hole_spacing);
        nodes(nodes_number(i) : a8roistis , 2) =  d/2 * cos(auto(1 : end - 1));
        nodes(nodes_number(i) : a8roistis , 3) =  d/2 * sin(auto(1 : end - 1));  
    end
    nodes(transpose(1 : find(nodes(:, 2), 1, 'last')), 4) = 1;       %id 1 = circle
   
    case 2
        perimeter = [0 0; 85 0; 85 50; 0 50];
         nodes(1,2:3) = 1;
         nodes(transpose(1 : find(nodes(:, 2), 1, 'last')), 4) = 1; 
         n_do_flat = x_nodes;
end

        
        





%% nodes generation perimeter

hold on

auto = dspace(perimeter(1, 1), perimeter(2, 1), n_do_flat, flat_bottom_spacing );
nodes(find(nodes(:,2), 1, 'last') + 1 : find(nodes(:,2), 1, 'last') + n_do_flat, 2) = auto(:);
nodes(find(nodes(:, 4), 1, 'last') + 1 : find(nodes(:, 4), 1, 'last') + n_do_flat, 4) = 2;       %id 2 = flat bottom

auto = dspace(perimeter(4, 1), perimeter(4, 2), n_r_flat, right_spacing);
nodes(find(nodes(:,2), 1, 'last') + 1 : find(nodes(:,2), 1, 'last') + n_r_flat, 2) = 85;
nodes(find(nodes(:,3), 1, 'last') + n_do_flat + 1 : find(nodes(:,3), 1, 'last') + n_do_flat + n_r_flat, 3) = auto(:);
nodes(find(nodes(:, 4), 1, 'last') + 1 : find(nodes(:, 4), 1, 'last') + n_r_flat, 4) = 3;       %id 3 = right

auto = dspace(perimeter(1, 1), perimeter(4, 2), n_le_flat, left_spacing );
nodes(find(nodes(:,2), 1, 'last') + 1 : find(nodes(:,2), 1, 'last') + n_le_flat, 2) = 0;
nodes(find(nodes(:,3), 1, 'last') + 1 : find(nodes(:,3), 1, 'last') + n_le_flat, 3) = auto(:); 
nodes(find(nodes(:, 4), 1, 'last') + 1 : find(nodes(:, 4), 1, 'last') + n_le_flat, 4) = 4;       %id 4 = left

auto = dspace(perimeter(2, 2), perimeter(2, 1), n_up_flat, flat_top_spacing );
nodes(find(nodes(:,2), 1, 'last') + 1 + n_le_flat : find(nodes(:,2), 1, 'last') + n_le_flat + n_up_flat, 2) = auto(:);
nodes(find(nodes(:,3), 1, 'last') + 1 : find(nodes(:,3), 1, 'last') + n_up_flat, 3) = 50;
nodes(find(nodes(:, 4), 1, 'last') + 1 : find(nodes(:, 4), 1, 'last') + n_up_flat, 4) = 5;        %id 5 = top

if epilogi == 2
    nodes(1, :) = [];
    nodes(:, 1) = nodes(:, 1) - 1;
end
    
%renumber and genaration of inner nodes


p = 1;
for i = 1 : size(nodes,1)
    if nodes(i, 4) == 4 
        renum_nodes(p, 1) = p;
        renum_nodes(p, 2 : 4) = nodes(i, 2 : 4); 
        p = p + 1;
    end
end
p = p -1;


% inner nodes


for i = 2 : x_nodes - 1

    random_x = transpose(dspace(nodes(i, 2), nodes(i + x_nodes + 2 * y_nodes, 2), y_nodes, x_spacing));
    random_y = transpose(dspace(nodes(i, 3), nodes(i + x_nodes + 2 * y_nodes, 3), y_nodes, y_spacing));
    
    %renumbering
    p = [p + 1 : p + size(random_x, 1)];
    renum_nodes(p(:), 1) = p(:);
    renum_nodes(p(:), 2) = random_x(:);
    renum_nodes(p(:), 3) = random_y(:);
    renum_nodes(p(1), 4) = nodes(i, 4);
    renum_nodes(p(end), 4) = nodes(i + x_nodes + 2 * y_nodes, 4); 
    p = p + size(random_x, 1)-1;      
end

thesi = p(1) + 1;
for i = 1 : size(nodes,1)
    if nodes(i, 4) ==  3
        renum_nodes(thesi, 1) = thesi;
        renum_nodes(thesi, 2 : 4) = nodes(i, 2 : 4); 
        thesi = thesi + 1;
    end
end

dupl_nodes = renum_nodes; 



%%
%mirror ws pros x

size_b4 = size(renum_nodes, 1);
k = 0;
%metra asfaleias

if epilogi == 2
    renum_nodes(1, 4) = 2;
end

for i = size_b4 + 1 : 2 * size_b4
        dupl_nodes(i , 1) =i - k;
        dupl_nodes(i , 2) = dupl_nodes(i - size_b4, 2);
        dupl_nodes(i , 3) = - dupl_nodes(i - size_b4, 3);
        dupl_nodes(i , 4) = dupl_nodes(i - size_b4, 4);
    
    if renum_nodes(i - size_b4, 4) == 2 || (renum_nodes(i- size_b4, 4) == 3 && (renum_nodes(i - size_b4, 3) == 0))
        k = k + 1;
        dupl_nodes(i, 1) = dupl_nodes(i - y_nodes * x_nodes, 1);
        
    else
        renum_nodes(i - k, 1) = renum_nodes(i - k - 1, 1) + 1;
        renum_nodes(i - k, 2) = renum_nodes(i - size_b4, 2);
        renum_nodes(i - k, 3) = - renum_nodes(i - size_b4, 3);
        renum_nodes(i - k, 4) = renum_nodes(i - size_b4, 4);
        
    end

end

l = k;
if epilogi == 2
    renum_nodes(1, 4) = 4;
end

%mirror ws pros y
size_b4 = size(renum_nodes, 1);
k = 0;
for i = size_b4 + 1 : 2 * size_b4
    
    if  renum_nodes(i - size_b4, 4) == 4 
        k = k + 1;
        
    else
        renum_nodes(i - k, 1) = renum_nodes(i - k - 1, 1) + 1;
        renum_nodes(i - k, 2) = - renum_nodes(i - size_b4, 2);
        renum_nodes(i - k, 3) = renum_nodes(i - size_b4, 3);
        renum_nodes(i - k, 4) = renum_nodes(i - size_b4, 4);
        
    end

end

size_b4 = size(dupl_nodes, 1);

for i = size_b4 + 1 : 2 * size_b4
 
    dupl_nodes(i , 1) = i - l;
    dupl_nodes(i , 2) = - dupl_nodes(i - size_b4, 2);
    dupl_nodes(i , 3) = dupl_nodes(i - size_b4, 3);
    dupl_nodes(i , 4) = dupl_nodes(i - size_b4, 4);
        
    if  dupl_nodes(i - size_b4, 4) == 4 
        dupl_nodes(i, 1) = dupl_nodes(i - 2 * y_nodes * x_nodes , 1);
        l = l + 1;
        
    end
    if ((dupl_nodes(i - size_b4, 4) == 2 && i > 3 * y_nodes * x_nodes)) 
        dupl_nodes(i, 1) = dupl_nodes(i - y_nodes * x_nodes, 1);
        l = l + 1;
        
        
    end
 end
dupl_nodes(end - y_nodes + 1, 1) = dupl_nodes(end - y_nodes + 1 - x_nodes * y_nodes, 1);   %lisi teleutaias stigmis
dupl_nodes(end - y_nodes + 2 : end, 1) = dupl_nodes(end - y_nodes + 2 : end, 1) - 1;       %same

g = 1;
for i = 1:size(dupl_nodes, 1)
    if dupl_nodes(i, 4) == 4
        pos(g) =  i;
        g = g+1;
    end
end




%fig = figure;  
scatter(dupl_nodes(:, 2), dupl_nodes(:, 3))
axis equal
savefig('Nodes plate.fig')
close(gcf)

%fig2 = figure;
scatter(dupl_nodes(:, 2), dupl_nodes(:, 3))
hold on


%%
%element creation

total_elements = ((x_nodes-1) * 2 * (y_nodes - 1)) * 4;
total_nodes = size(renum_nodes,1);

elements = struct('id', [], 'Nodes', [], 'CG', [], 't', [], 'A', []);
elements.Nodes = zeros(total_elements, 3);
%kappa = 1;
position = 1;  % metritis gia ton ari8mo twn element / epanalipsi (2)
active = 1;     % metritis gia tin 8esi sta nodes
tessera = 1;
for i = 1 : total_elements / 2
    sit1 = sqrt((dupl_nodes(active + y_nodes + 1, 2) - dupl_nodes(active, 2))^ 2 +(dupl_nodes(active + y_nodes + 1, 3) - dupl_nodes(active, 3))^ 2);
    sit2 = sqrt((dupl_nodes(active + 1, 2) - dupl_nodes(active + y_nodes, 2))^ 2 +(dupl_nodes(active + 1, 3) - dupl_nodes(active + y_nodes, 3))^ 2);
    if sit1 >= sit2
        if i > total_elements / 8 && i <=3 * total_elements / 8
            elements.id(position) = position;
            elements.Nodes(position, 1) = dupl_nodes(active, 1);
            elements.Nodes(position, 3) = dupl_nodes(active + y_nodes, 1);
            elements.Nodes(position, 2) = dupl_nodes(active + y_nodes + 1, 1);
            position = position + 1;
        
            elements.id(position) = position;
            elements.Nodes(position, 1) = dupl_nodes(active, 1);
            elements.Nodes(position, 3) = dupl_nodes(active + y_nodes + 1, 1);
            elements.Nodes(position, 2) = dupl_nodes(active + 1, 1);
            position = position + 1;
        else
            
            elements.id(position) = position;
            elements.Nodes(position, 1) = dupl_nodes(active, 1);
            elements.Nodes(position, 2) = dupl_nodes(active + y_nodes, 1);
            elements.Nodes(position, 3) = dupl_nodes(active + y_nodes + 1, 1);
            position = position + 1;
        
            elements.id(position) = position;
            elements.Nodes(position, 1) = dupl_nodes(active, 1);
            elements.Nodes(position, 2) = dupl_nodes(active + y_nodes + 1, 1);
            elements.Nodes(position, 3) = dupl_nodes(active + 1, 1);
            position = position + 1;
        end
    else
        if i > total_elements / 8 && i <=3 * total_elements / 8
            elements.id(position) = position;
            elements.Nodes(position, 1) = dupl_nodes(active, 1);
            elements.Nodes(position, 3) = dupl_nodes(active + y_nodes, 1);
            elements.Nodes(position, 2) = dupl_nodes(active + 1, 1);
            position = position + 1;
        
            elements.id(position) = position;
            elements.Nodes(position, 1) = dupl_nodes(active + y_nodes, 1);
            elements.Nodes(position, 3) = dupl_nodes(active + y_nodes + 1, 1);
            elements.Nodes(position, 2) = dupl_nodes(active + 1, 1);
            position = position + 1;
        else
            
            elements.id(position) = position;
            elements.Nodes(position, 1) = dupl_nodes(active, 1);
            elements.Nodes(position, 2) = dupl_nodes(active + y_nodes, 1);
            elements.Nodes(position, 3) = dupl_nodes(active + 1, 1);
            position = position + 1;
        
            elements.id(position) = position;
            elements.Nodes(position, 1) = dupl_nodes(active + y_nodes, 1);
            elements.Nodes(position, 2) = dupl_nodes(active + y_nodes + 1, 1);
            elements.Nodes(position, 3) = dupl_nodes(active + 1, 1);
            position = position + 1;
        end
        
    end
    
    active = active + 1;
    
    if mod(i,(total_elements/8)) == 0
        active = active + y_nodes;
    end
    if mod(active,y_nodes) == 0
        active = active + 1;
    end
   line(renum_nodes(elements.Nodes(position - 2, [1, 2, 3, 1]), 2), renum_nodes(elements.Nodes(position - 2, [1, 2, 3, 1]), 3));
   line(renum_nodes(elements.Nodes(position - 1, [1, 2, 3, 1]), 2), renum_nodes(elements.Nodes(position - 1, [1, 2, 3, 1]), 3));


end

elements.id = transpose(elements.id);

axis equal
title('Meshed Plate')
xlabel('mm')
ylabel('mm')
savefig('meshed_plaka.fig')
close(gcf)
%%
% Final Nodes Definition
Nodes = struct('id', [], 'coords', [], 'BCs', [], 'F', [], 'Flag', []);
Nodes.id = renum_nodes(:, 1);
Nodes.coords = renum_nodes(:, 2 : 3);
Nodes.F = zeros(size(Nodes.coords, 1), 2);
Nodes.BCs = zeros(size(Nodes.coords, 1), 2);
Nodes.Flag = renum_nodes(:, 4);


%%
% Element Rest attributes
elements.CG = zeros(size(elements.Nodes, 1), 2); elements.t = zeros(size(elements.Nodes, 1), 1); elements.A = zeros(size(elements.Nodes, 1), 1);

for i = 1 : size(elements.Nodes, 1)
   elements.CG(i,1) = (Nodes.coords(elements.Nodes(i, 1), 1) + Nodes.coords(elements.Nodes(i, 2), 1) + Nodes.coords(elements.Nodes(i, 3), 1)) / 3;
   elements.CG(i,2) = (Nodes.coords(elements.Nodes(i, 1), 2) + Nodes.coords(elements.Nodes(i, 2), 2) + Nodes.coords(elements.Nodes(i, 3), 2)) / 3;
   elements.t(i, 1) = t;
   elements.A(i, 1) = 1 / 2 * det([Nodes.coords(elements.Nodes(i, 1), 1) Nodes.coords(elements.Nodes(i, 1), 2) 1;
                            Nodes.coords(elements.Nodes(i, 2), 1) Nodes.coords(elements.Nodes(i, 2), 2) 1;
                            Nodes.coords(elements.Nodes(i, 3), 1) Nodes.coords(elements.Nodes(i, 3), 2) 1]);
end

%scatter(elements.CG(:, 1), elements.CG(:, 2), "filled", "g")


%%
%Forces (Dark side)
q = 1;
while q ~= 0
    stress_type = listdlg('PromptString', 'Define Load case', 'SelectionMode', 'single', 'ListString', {'Tensile', 'Bend'});
    if isnumeric(stress_type) == 1
        q = 0;
    else
        waitfor(msgbox('A selection is needed'));
    end
end


if stress_type == 1
    q = 0;
    while q == 0
        Ften = inputdlg('Define Tensile Force measure (F)', 'Tensile Force input', dims);
        Ften = str2double(Ften);
        F_nom = Ften;
        if isnumeric(Ften) == 1 && Ften >=0
            q = 1;
        else
            waitfor(msgbox('wrong input, Please give a positivve number'))
        end
    end
    [Nodes.F, answer] = Force_Selection2D(Nodes.id, Nodes.coords, stress_type, Nodes.Flag, Nodes.F, Ften, leny, y_nodes);
else
    q = 0;
    while q == 0
        Mben = inputdlg('Define Bending Torque measure (Nmm)', 'Bending Torque input', dims);
        Mben = str2double(Mben);
        F_nom = Mben;
        if isnumeric(Mben) == 1 && Mben >=0
            q = 1;
        else
            waitfor(msgbox('wrong input, Please give a positive number'))
        end
    end
    [Nodes.F, answer] = Force_Selection2D(Nodes.id, Nodes.coords, stress_type, Nodes.Flag, Nodes.F, Mben, leny, y_nodes, t);
end

%% BCs
Nodes.BCs = Boundary_Condition_Selection2D(Nodes.coords, Nodes.BCs, stress_type, answer);
%% Final geometry
hold on
patch('Faces', elements.Nodes, 'Vertices', Nodes.coords, 'FaceColor', 'none')
scatter(Nodes.coords(find(Nodes.BCs(:, 1) == 1), 1), Nodes.coords(find(Nodes.BCs(:, 1) == 1), 2), 'filled', 'green')
scatter(Nodes.coords(find(Nodes.BCs(:, 2) == 1), 1), Nodes.coords(find(Nodes.BCs(:, 2) == 1), 2), 'filled', 'magenta')
scatter(Nodes.coords(find(Nodes.F(:, 1) ~= 0), 1), Nodes.coords(find(Nodes.F(:, 1) ~= 0), 2), 'filled', 'yellow')
quiver(Nodes.coords(find(Nodes.F(:, 1) ~= 0), 1), Nodes.coords(find(Nodes.F(:, 1) ~= 0), 2), Nodes.F(find(Nodes.F(:, 1) ~= 0), 1), Nodes.F(find(Nodes.F(:, 1) ~= 0), 2))

axis equal
hold off


%% TXT Files

%create txt file for nodes
fID = fopen('AEM6101_Nodes_Plate.txt', 'w');
fprintf(fID,'%3s,%13.5s,%13.5s,%5s,%5s,%8.8s,%8.8s,%3s \n', 'id', 'x', 'y', 'BCx', 'BCy', 'Fx', 'Fy', 'Flags');
for i = 1 : size(Nodes.id)
    fprintf(fID,'%3d,%13.5d,%13.5d,%5d,%5d,%8.8d,%8.8d,%3d \n', Nodes.id(i), Nodes.coords(i,1), Nodes.coords(i, 2),Nodes.BCs(i, 1), Nodes.BCs(i, 2), Nodes.F(i, 1), Nodes.F(i, 2),Nodes.Flag(i));
end
fclose(fID);

%create txt file for elements
fID = fopen('AEM6101_Elements_Plate.txt', 'w');
fprintf(fID,'%5s, %4s, %4s, %4s, %8.2s, %8.3s, %8.3s, %4s \n', 'id', 'N1', 'N2', 'N3', 'A', 'CGx', 'CGy', 't');
for i = 1 : size(elements.id)
    fprintf(fID,'%5d, %4d, %4d, %4d, %8.4d, %8.2d, %8.2d, %4d \n', elements.id(i, 1), elements.Nodes(i, 1), elements.Nodes(i, 2), elements.Nodes(i, 3), elements.A(i, 1), elements.CG(i, 1), elements.CG(i, 2), elements.t(i, 1));
end
fclose(fID);

%create txt file for data
fID = fopen('AEM6101_data_Plate.txt', 'w');
fprintf(fID,'%8.2s, %2.2s, %8.3s, %5.2s, %5.2s, %5.2s, %5.2s, %5.2s, %5.2s, %5.2s, %5.2s \n', 'E', 'v', 'G', 'lenx', 'leny', 'x_nodes', 'y_nodes', 'Force_type', 'Fnom', 'd', 't');
fprintf(fID,'%8.2d, %2.2d, %8.3d, %5.2d, %5.2d, %5.2d, %5.2d, %5.2d, %5.2d, %5.2d, %5.2d \n', E, v, G, 2 * lenx, 2 * leny, x_nodes, y_nodes, stress_type, F_nom, d, t);
fclose(fID);



%% Final Options

choice = questdlg('Are you satisfied with your mesh?');
if choice == "Yes"
    choice2 = questdlg('Iniciate solver?');
    if choice2 == "Yes"
        AEM6101_SOLVER
    else
        msgbox('Pre-Processor terminated')
        return
    end
else
    choice3 = questdlg('Start meshing over?');
    if choice3 == "Yes"
        AEM6101_PRE_PROCESSOR
    else
        msgbox('Pre-Processor terminated')
        return
    end
end
        

