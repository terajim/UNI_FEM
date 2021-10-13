%geranos solver
%Terzis Dimitrios AEM6101

%Check if nodes file exist

if exist('./AEM6101_nodes_ger_truss.txt','file')==2
     nodes_ger=importdata('AEM6101_nodes_ger_truss.txt',',');
else
     waitfor(msgbox('Nodes file not found!'))
     error('Nodes file not found')
end

%check if elements file exist
if exist('./AEM6101_elements_ger_truss.txt','file')==2
     elements_ger=importdata('AEM6101_elements_ger_truss.txt',',');
else
     waitfor(msgbox('Elements file not found!'))
     error('Elements file not found')
end

%check if force file exist
if exist('./AEM6101_force_ger_truss.txt','file')==2
     force_ger=importdata('AEM6101_force_ger_truss.txt',',');
else
     waitfor(msgbox('Forces file not found!'))
     error('Forces file not found')
end

%check if bc file exist
if exist('./AEM6101_BC_ger_truss.txt','file')==2
     bc_nodes_ger=importdata('AEM6101_BC_ger_truss.txt',',');
else
     waitfor(msgbox('Boundary Conditions file not found!'))
     error('Boundary conditions file not found')
end

figure
str=ones(size(nodes_ger,1),1);
    for j=1:size(nodes_ger,1)
        str(j,1)=j;
        scatter3(nodes_ger(j,2),nodes_ger(j,3),nodes_ger(j,4),'red','filled')
        text(nodes_ger(j,2),nodes_ger(j,3),nodes_ger(j,4),{'node',num2str(str(j,1))})
        hold on
    end
    title('Geranos') %title of plots
    xlabel('x-axis [mm]')   %specify units of measurement
    ylabel('y-axis [mm]')   
    zlabel('z-axis [mm]')
    for i=1:8*Storeys+18
        line(nodes_ger(elements(i,2:3),2),nodes_ger(elements(i,2:3),3),nodes_ger(elements(i,2:3),4),'LineWidth',3)
    end 
    for i=8*Storeys+19:size(elements_ger,1)-3
        line(nodes_ger(elements(i,2:3),2),nodes_ger(elements(i,2:3),3),nodes_ger(elements(i,2:3),4),'LineWidth',1)
    end 
    for i=size(elements,1)-3:size(elements_ger,1)
        line(nodes_ger(elements_ger(i,2:3),2),nodes_ger(elements_ger(i,2:3),3),nodes_ger(elements_ger(i,2:3),4),'Color','b','LineWidth',1)
    end    
    axis equal
    camproj('orthographic')
    view([-8 -15 1])

%initialize K stiffnes matrix of cranecu
K_stiff=zeros(size(nodes_ger,1)*3); %number of nodes_ger * degrees of freedom
K_Stiff_global=zeros(size( nodes_ger,1)*3,size( nodes_ger,1)*3);
elements_length=zeros(size(elements_ger,1));
for i=1:size( elements_ger,1)
    elements_length(i,1)=norm( nodes_ger( elements_ger(i,3),2:4)- nodes_ger( elements_ger(i,2),2:4));
    %local coordinate system 
    l=( nodes_ger( elements_ger(i,3),2)- nodes_ger( elements_ger(i,2),2))/elements_length(i,1);
    m=( nodes_ger( elements_ger(i,3),3)- nodes_ger( elements_ger(i,2),3))/elements_length(i,1);
    n=( nodes_ger( elements_ger(i,3),4)- nodes_ger( elements_ger(i,2),4))/ elements_length(i,1);
    ke=(( elements_ger(i,4)* elements_ger(i,5))/elements_length(i,1))*[1 -1;-1 1];
    T=[l m n 0 0 0; 0 0 0 l m n]; %matrix of directions
    Ke=T'*ke*T;
    
    if elements_ger(i,1)>elements_ger(i,2)
         a=3*elements_ger(i,2)-2;
         z=3*elements_ger(i,3)-2;
    else
         a=3*elements_ger(i,2)-2;
         z=3*elements_ger(i,3)-2;
    end
    
    K_Stiff_global(a:a+2,a:a+2)=K_Stiff_global(a:a+2,a:a+2)+Ke(1:3,1:3);
    K_Stiff_global(a:a+2,z:z+2)=K_Stiff_global(a:a+2,z:z+2)+Ke(1:3,4:6);
    K_Stiff_global(z:z+2,a:a+2)=K_Stiff_global(z:z+2,a:a+2)+Ke(4:6,1:3);
    K_Stiff_global(z:z+2,z:z+2)=K_Stiff_global(z:z+2,z:z+2)+Ke(4:6,4:6);    
end         

%penalty method
K=K_Stiff_global;
for i=1:size(bc_nodes_ger,1)
    if bc_nodes_ger(i,2)==1
        K(3*i-2,3*i-2)=K(3*i-2,3*i-2)+10^8;
    end
    if bc_nodes_ger(i,3)==1
        K(3*i-1,3*i-1)=K(3*i-1,3*i-1)+10^8;
    end
    if bc_nodes_ger(i,4)==1
        K(3*i,3*i)=K(3*i,3*i)+10^8;
    end
end


%kanw ton pinaka f, mia stilli
F=reshape(force_ger(:,2:4)',[],1); 

%check if K is invertible
D=det(K);
if D==0
    msgbox('Error! Stiffness Matrix not invertible! Check Boundary Conditions')
    error('Error! K not invertible')
end


%Solve F=K*u
u=inv(K)*F; %nodes displacements
%Calculate Reaction Forces
Rs=K_Stiff_global*u; 
U=vec2mat(u,3);
strain=zeros(size(elements_ger,1),1);
stress=zeros(size(elements_ger,1),1);
nodes_newcord=zeros(size(nodes_ger,1),4);
for i=1:size(nodes_ger,1)
    nodes_newcord(i,1)=nodes_ger(i,1);
    nodes_newcord(i,2)=nodes_ger(i,2)+U(i,1);
    nodes_newcord(i,3)=nodes_ger(i,3)+U(i,2);
    nodes_newcord(i,4)=nodes_ger(i,4)+U(i,3);
end

for i=1:size(elements_ger,1)
    elements_newlength(i,1)=norm(nodes_newcord(elements_ger(i,3),2:4)-nodes_newcord(elements_ger(i,2),2:4));
    %calculate strains of each element
    strain(i,1)=(elements_newlength(i,1)-elements_length(i,1))/(elements_length(i,1));
    stress(i,1)=210000*strain(i,1);
end


%convert Reactions to Matrix for 3dof
Reactions=vec2mat(Rs,3);


%calculate total reaction on each node
Rtotal(:,1)=sqrt(Reactions(:,1).^2+Reactions(:,2).^2+Reactions(:,3).^2);
Reactions=[Reactions,Rtotal]; 



%create txt file for nodes first coordinates
fID=fopen('AEM6101_Nodes_ger.txt','w');
fprintf(fID,'%3s,%6s,%6s,%6s \r\n','id','x','y','z');
for i=1:size(nodes_ger,1)
        fprintf(fID,'%3d,%6d,%6d,%6d \r\n',nodes_ger(i,1),nodes_ger(i,2),nodes_ger(i,3),nodes_ger(i,4));
end
fclose(fID);

%create txt file for nodes new coordinates
fID=fopen('AEM6101_NodesNewCord_ger.txt','w');
fprintf(fID,'%3s,%15s,%15s,%15s \r\n','id','newx','newy','newz');
for i=1:size(nodes_newcord,1)
    fprintf(fID,'%3d,%15d,%15d,%15d \r\n',nodes_newcord(i,1),nodes_newcord(i,2),nodes_newcord(i,3),nodes_newcord(i,4));
end
fclose(fID);


%create txt file for boundary conditions
fID=fopen('AEM6101_BC_ger.txt','w');
fprintf(fID,'%3s,%5s,%5s,%5s \r\n','id','bc_x','bc_y','bc_z');
for i=1:size(bc_nodes_ger,1)
    fprintf(fID,'%3d,%5d,%5d,%5d \r\n',bc_nodes_ger(i,1),bc_nodes_ger(i,2),bc_nodes_ger(i,3),bc_nodes_ger(i,4));
end
fclose(fID);

%create txt file for forces
fID=fopen('AEM6101_Forces_ger.txt','w');
fprintf(fID,'%3s,%5s,%5s,%5s \r\n','id','Fx','Fy','Fz');
for i=1:size(force_ger,1)
    fprintf(fID,'%3d,%5d,%5d,%5d \r\n',force_ger(i,1),force_ger(i,2),force_ger(i,3),force_ger(i,4));
end
fclose(fID);

%create txt file for elements
fID=fopen('AEM6101_Elements_ger.txt','w');
fprintf(fID,'%5s,%4s,%4s,%8s \r\n','id','Node1','Node2','A');
for i=1:size(elements_ger,1)
    fprintf(fID,'%5d,%4d,%4d,%8.4d,%4d \r\n',elements_ger(i,1),elements_ger(i,2),elements_ger(i,3),elements_ger(i,4),elements_ger(i,5));
end
fclose(fID);


%create txt file for elements length
fID=fopen('AEM6101_Elementslength_ger.txt','w');
fprintf(fID,'%15s \r\n','length');
for i=1:size(elements_length,1)
    fprintf(fID,'%15d \r\n',elements_length(i,1));
end
fclose(fID);

%create txt file for elements new length
fID=fopen('AEM6101_ElementsNewlength_ger.txt','w');
fprintf(fID,'%15s \r\n','newlength');
for i=1:size(elements_newlength,1)
    fprintf(fID,'%15d \r\n',elements_newlength(i,1));
end
fclose(fID);


%create txt file for Strains 
fID=fopen('AEM6101_Strains_ger.txt','w');
fprintf(fID,'%12.7s \r\n','strains');
for i=1:size(strain,1)
    fprintf(fID,'%12.4d \r\n',strain(i,1));
end
fclose(fID);

%create txt file for Displacements
fID=fopen('AEM6101_Displacement_ger.txt','w');
fprintf(fID,'%11.4s,%11.4s,%11.4s \r\n','Ux','Uy','Uz');
for i=1:size(U,1)
    fprintf(fID,'%11.4d,%11.4d,%11.4d \r\n',U(i,1),U(i,2),U(i,3));
end
fclose(fID);


%create Reactions txt file
fID=fopen('AEM6101_Reactions_ger.txt','w');
fprintf(fID,'%15.2s,%15.2s,%15.2s,%15.2s \r\n','Rx','Ry','Rz','Rtotal');
for i=1:size(Reactions,1)
    fprintf(fID,'%15.2d,%15.2d,%15.2d,%15.2d \r\n',Reactions(i,1),Reactions(i,2),Reactions(i,3),Reactions(i,4));
end
fclose(fID);


%create Stresses txt file
fID=fopen('AEM6101_Stresses_ger.txt','w');
fprintf(fID,'%12.8s \n','stresses');
for i=1:size(stress,1)
    fprintf(fID,'%12.4d \n',stress(i,1));
end
fclose(fID);
 


str=ones(8,1);
    for j=1:size(nodes_newcord,1)
        str(j,1)=j;
        scatter3(nodes_newcord(j,2),nodes_newcord(j,3),nodes_newcord(j,4),'red','filled')
        text(nodes_newcord(j,2),nodes_newcord(j,3),nodes_newcord(j,4),{'node',num2str(str(j,1))})
        hold on
    end

    for i=1:size(elements_ger,1)
        line(nodes_newcord(elements_ger(i,2:3),2),nodes_newcord(elements_ger(i,2:3),3),nodes_newcord(elements_ger(i,2:3),4),'LineWidth',1)
    end

    
    msgbox('Congratulations! Sucess! Good job! Proceed to Post processor')
    