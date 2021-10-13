function[BC_NODES] = Boundary_Condition_Selection2D(nodes, bc_nodes, stress_type, answer)


% ask user if wanna change boundary Contitions
 
                %%%%-----------geometry pop up--------------%%%%
    fig = figure;
    scatter(nodes(:, 1), nodes(:, 2))
    axis equal
    hold on
                        %%%%-------------------------%%%%


if answer == 2  
    q1 = 1;
    while q1 == 1
         list = listdlg('PromptString','Please select a Boundary Condition Style','SelectionMode','single','ListString',{'x,y','x','y'});
         k = 1;
         while k~= 0
             if isnumeric(list) == 1
                 k = 0;
             else
                 waitfor(msgbox('A selection is needed'));
             end
         end
         hold on
         if list == 1                 
             hold on
             shg
             dcm_obj = datacursormode(fig);
             set(dcm_obj,'DisplayStyle','window',...
             'SnapToDataVertex', 'on', 'Enable', 'on')
             waitfor(msgbox('Please select a node to lock its x, y'))
             waitforbuttonpress
             c_info = getCursorInfo(dcm_obj);
             for i=1:size(nodes,1)
                if sum(c_info.Position == nodes(i, :)) == 2
                    bc_nodes(i, :) = 1;
                    scatter(c_info.Position(1), c_info.Position(2), 'filled', 'm')
                end
             end
             datacursormode off
             q1=listdlg('PromptString','Do you wish to add more?','SelectionMode','single','ListString',{'Yes','No'}); 
             k = 1;
             while k~= 0
                if isnumeric(q1) == 1
                    k = 0;
                else
                    waitfor(msgbox('A selection is needed'));
                end
             end
         elseif list == 2
             hold on
             shg
             dcm_obj = datacursormode(fig);
             set(dcm_obj,'DisplayStyle','window',...
             'SnapToDataVertex','on','Enable','on')
             waitfor(msgbox('Please select a node to lock its x'))
             waitforbuttonpress
             c_info = getCursorInfo(dcm_obj);
             for i = 1:size(nodes,1)
                if sum(c_info.Position == nodes(i, :)) == 2
                    bc_nodes(i, 1) = 1;
                    scatter(c_info.Position(1), c_info.Position(2), 'filled', 'b')
                end
             end
             datacursormode off
             q1=listdlg('PromptString','Do you wish to add more?','SelectionMode','single','ListString',{'Yes','No'}); 
             k = 1;
             while k~= 0
                if isnumeric(q1) == 1
                    k = 0;
                else
                    waitfor(msgbox('A selection is needed'));
                end
            end
         else 
             hold on
             shg
             dcm_obj = datacursormode(fig);
             set(dcm_obj,'DisplayStyle','window',...
             'SnapToDataVertex','on','Enable','on')
             waitfor(msgbox('Please select a node to lock its y'))
             waitforbuttonpress
             c_info = getCursorInfo(dcm_obj);
             for i=1:size(nodes,1)
                if sum(c_info.Position == nodes(i, :)) == 3
                    bc_nodes(i, 2)=1;
                    scatter(c_info.Position(1), c_info.Position(2), 'filled', 'r')
                end
             end
             datacursormode off
             q1=listdlg('PromptString','Do you wish to add more?','SelectionMode','single','ListString',{'Yes','No'});
             while k~= 0
                if isnumeric(q1) == 1
                    k = 0;
                else
                    waitfor(msgbox('A selection is needed'));
                end
            end
         end
         
    end
    
else
    
    for i=1:size(nodes, 1)
        if nodes(i, 1) == - 85
            hold on
            if  nodes(i, 2) == 0
                bc_nodes(i, 2) = 1;
                scatter(nodes(i, 1), nodes(i, 2), 'filled', 'r') 
            end
            hold on
            bc_nodes(i, 1) = 1;
            scatter(nodes(i, 1), nodes(i, 2), 'filled', 'g')
        end
    end
end
hold off

%scatter(nodes(bc_nodes == 1, 1), nodes(bc_nodes == 1, 1), 'filled', 'r')
 
 hold off
 BC_NODES = bc_nodes;
 close all
 
