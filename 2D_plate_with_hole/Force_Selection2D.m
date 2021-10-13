function[FORCE,answer] = Force_Selection2D(id, nodes, stress_type, flags, force, Force_meas, leny, y_nodes, t)

q1 = 1;
% ask user if wanna change boundary Contitions
k = 1;
while k ~= 0
    answer = listdlg('PromptString', 'Use Default force and BCs?', 'SelectionMode', 'single', 'ListString', {'Yes', 'No'});
    if isnumeric(answer) == 1
        k = 0;
    else
        waitfor(msgbox('A selection is needed'))
    end
end 
                %%%%-----------geometry pop up--------------%%%%
    fig = figure;
    scatter(nodes(:, 1), nodes(:, 2))
    axis equal
    hold on
    %line(nodes(elements(:, [1, 2, 3, 1]), 2), nodes(elements(:, [1, 2, 3, 1]), 3));
    %axis equal-
    
if answer == 2
    
    while q1 == 1
        k = 1;
        while k~= 0
            hold on
            shg
            dcm_obj = datacursormode(fig);
            set(dcm_obj,'DisplayStyle','window',...
            'SnapToDataVertex', 'on', 'Enable', 'on')
            waitfor(msgbox('Please select an edge node to apply force'));
            waitforbuttonpress
            c_info = getCursorInfo(dcm_obj);
            for i=1:size(nodes,1) 
                if sum(c_info.Position == nodes(i, :)) == 2
                    flag = i; 
                end                
            end
            if flags(flag) == 0 || flags(flag) == 1 || flags(flag) == 2 || flags(flag) == 4 || flags(flag) == 5
                waitfor(msgbox('Wrong input. Please select an edge node.'));
            else
                k = 0;
            end
        end
        multiple = listdlg('PromptString', 'Do you wish to select the whole column at once?', 'SelectionMode', 'single', 'ListString', {'Yes', 'No'});
        if multiple == 1
                NodesWithForce = id(flags == flags(flag));
                if flag < (size(nodes, 1) / 2) 
                    NodesWithForce = NodesWithForce(1 : 2 * y_nodes - 1);
                else
                    NodesWithForce = NodesWithForce(2 * y_nodes  : end);
                end
                NodesWithForce(:, 2:3) = nodes(NodesWithForce(:), :);
                NodesWithForce(:, 3) = sort(NodesWithForce(:, 3));
                m = NodesWithForce(1: y_nodes , 1);
                NodesWithForce(1: y_nodes - 1, 1) = sort(NodesWithForce(y_nodes + 1 : end, 1), 'descend');
                NodesWithForce(y_nodes : end, 1) = m;
                scatter(NodesWithForce(:, 2),NodesWithForce(:, 3), 'filled', 'g')
                
        else
            
            force(i, :) = 1;
            scatter(c_info.Position(1), c_info.Position(2), 'filled', 'g')
            waitfor(msgbox('Be carefull at the BC selection, defaults will not work for a single node tension'));
        
        end
        datacursormode off
        k = 1;
        while k~= 0
            q1 = listdlg('PromptString','Do you wish to add more?','SelectionMode','single','ListString',{'Yes','No'}); 
            if isnumeric(q1) == 1
                k = 0;
            else
                waitfor(msgbox('A selection is needed'));
            end
        end
    end
    
else
    
    hold on
    NodesWithForce = id(flags == 3);
    NodesWithForce = NodesWithForce(1 : 2 * y_nodes - 1);
    NodesWithForce(:, 2:3) = nodes(NodesWithForce(:), :);
    NodesWithForce(:, 3) = sort(NodesWithForce(:, 3));
    flag = NodesWithForce(1: y_nodes , 1);
    NodesWithForce(1: y_nodes - 1, 1) = sort(NodesWithForce(y_nodes + 1 : end, 1), 'descend');
    NodesWithForce(y_nodes : end, 1) = flag;
    scatter(NodesWithForce(:, 2),NodesWithForce(:, 3), 'filled', 'g')
    
    %NodesWithForce
    
end
Hel = NodesWithForce(2:end, 3) - NodesWithForce(1 :end-1, 3); % Heights of Elements on side that is loaded
if stress_type == 1                        %1 = tesnile
    Fel = Force_meas * Hel ./ sum(Hel);    % Total Force on each Element FOR TENSILE 
else
    Del = NodesWithForce(:, 3) - ones(size(NodesWithForce, 1), 1) * NodesWithForce(ceil((2 * y_nodes - 1) / 2), 3); %Distance of each node from centerline
    Sel = Force_meas .* Del  / (t *(2 * leny) ^3 / 12);      % Stress on each Node. apo analitiki lisi
    Fel = (Sel(2 : end) + Sel(1 : end - 1)).* Hel(:)./ 2;  % Total Force on each Element BENDING

end


% ForcesRS: (:,2)Array. left column is x dir, right is ydir. It is the Force
% in each Node on the loaded side.
ForcesRS(1,:) = [Fel(1, :) / 2, 0];
ForcesRS(2 : 2* y_nodes - 2, :) = [(Fel(2 : end, 1) + Fel(1 : end - 1, 1))/2, zeros(2* y_nodes - 3, 1)];
ForcesRS(2 * y_nodes - 1, :) = [Fel(end, :)/2, 0];
if stress_type == 2
    sum(ForcesRS(:, 1).* Del(:))
    ForcesRS = ForcesRS*4.0128;
    sum(ForcesRS(:, 1).* Del(:))
end
 
for i = 1: size(force, 1)
    for j = 1: size(NodesWithForce, 1)
        if id(i) == NodesWithForce(j, 1)
            force(i,:) = ForcesRS(j, :);
        end
    end
end
FORCE = force;
close all
end

 


 
 