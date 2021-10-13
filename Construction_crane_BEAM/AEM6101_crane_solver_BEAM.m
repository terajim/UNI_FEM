%Solver Beam Geranos
%Terzis Dimitrios AEM6101

%Check if nodes file exist

ID=listdlg('PromptString','Select solver purpose','SelectionMode','single','ListString',{'Beam cube','Beam crane'});
if ID==2

    if exist('./AEM6101_nodes_ger_beam.txt','file')==2
        nodes_beam=importdata('AEM6101_nodes_ger_beam.txt',',');
    else
        waitfor(msgbox('Nodes file not found!'))
        error('Nodes file not found')
    end


%check if elements file exist
    if exist('./AEM6101_elements_ger_beam.txt','file')==2
        elements_beam=importdata('AEM6101_elements_ger_beam.txt',',');
    else
        waitfor(msgbox('Elements file not found!'))
        error('Elements file not found')
    end

%check if force file exist
    if exist('./AEM6101_force_ger_beam.txt','file')==2
         force_beam=importdata('AEM6101_force_ger_beam.txt',',');
    else
        waitfor(msgbox('Forces file not found!'))
        error('Forces file not found')
    end


%check if bc file exist
    if exist('./AEM6101_BC_ger_beam.txt','file')==2
        bc_nodes_beam=importdata('AEM6101_BC_ger_beam.txt',',');
    else
        waitfor(msgbox('Boundary Conditions file not found!'))
        error('Boundary conditions file not found')
    end
    
    
elseif ID==1
    
    
    
    if exist('./AEM6101_cube_nodes_beam.txt','file')==2
        nodes_beam=importdata('AEM6101_cube_nodes_beam.txt',',');
    else
        waitfor(msgbox('Nodes file not found!'))
        error('Nodes file not found')
    end


%check if elements file exist
    if exist('./AEM6101_cube_elements_beam.txt','file')==2
        elements_beam=importdata('AEM6101_cube_elements_beam.txt',',');
    else
        waitfor(msgbox('Elements file not found!'))
        error('Elements file not found')
    end

%check if force file exist
    if exist('./AEM6101_cube_force_beam.txt','file')==2
         force_beam=importdata('AEM6101_cube_force_beam.txt',',');
    else
        waitfor(msgbox('Forces file not found!'))
        error('Forces file not found')
    end


%check if bc file exist
    if exist('./AEM6101_BC_cube_beam.txt','file')==2
        bc_nodes_beam=importdata('AEM6101_BC_cube_beam.txt',',');
    else
        waitfor(msgbox('Boundary Conditions file not found!'))
        error('Boundary conditions file not found')
    end
    
else
    error('Wrong input,please start over!')
end


v=0.3; %poisson's ratio

G=210000/(2*(1+v)); %shear modulus
E=210000;
%initialize K stiffnes matrix of cube , 1 for beam 2 for truss
K_stiff1=zeros((size(nodes_beam,1)-1)*6); %BEAMS %number of nodes * degrees of freedom
K_stiff2=zeros((size(nodes_beam,1)-1)*6); %TRUSS



for i=1:size(elements_beam,1)
    if i==59
        ke=zeros(12); %number of nodes(4) * dof(3) , for elements 59-62 (TRUSS)
    end
    if i<59
         X1=nodes_beam(elements_beam(i,2),2); %x 1st node, element i
         Y1=nodes_beam(elements_beam(i,2),3); %y 1st node, element i
         Z1=nodes_beam(elements_beam(i,2),4); %z 1st node, element i
         X2=nodes_beam(elements_beam(i,3),2); %x 2nd node, element i
         Y2=nodes_beam(elements_beam(i,3),3); %y 2nd node, element i
         Z2=nodes_beam(elements_beam(i,3),4); %z 2nd node, element i
         X3=nodes_beam(elements_beam(i,4),2); %x node 33 (0)
         Y3=nodes_beam(elements_beam(i,4),3); %y node 33 (0)
         Z3=nodes_beam(elements_beam(i,4),4); %z node 33 (~=0) 
   
         elements_beam_length(i,1)=norm(nodes_beam(elements_beam(i,3),:)-nodes_beam(elements_beam(i,2),:)); % distance between: (2nd node's x of element i)-(1st node's x of element i)
         
         A123=sqrt(((Y2-Y1)*(Z3-Z1)-(Y3-Y1)*(Z2-Z1))^2+((Z2-Z1)*(X3-X1)-(Z3-Z1)*(X2-X1))^2+((X2-X1)*(Y3-Y1)-(X3-X1)*(Y2-Y1))^2);
         
         R=sqrt(elements_beam(i,5)/pi); %beam's radi
         a=elements_beam_length(i,1)/2;
         Iy=(pi*R^4)/4;
         Iz=Iy;
         J=2*Iy;   %polar moment of inertia
         
         %cos, sin etc for directions x y z
         lx=(X2-X1)/elements_beam_length(i,1); 
         mx=(Y2-Y1)/elements_beam_length(i,1); 
         nx=(Z2-Z1)/elements_beam_length(i,1); 
         lz=((Y2-Y1)*(Z3-Z1)-(Y3-Y1)*(Z2-Z1))/A123;
         mz=((Z2-Z1)*(X3-X1)-(Z3-Z1)*(X2-X1))/A123;
         nz=((X2-X1)*(Y3-Y1)-(X3-X1)*(Y2-Y1))/A123;
         ly=mz*nx-nz*mx;
         my=nz*lx-lz*nx;
         ny=lz*mx-mz*lx;
         
         %direction matrix T3 
         T3=[lx mx nx;ly my ny;lz mz nz];
         % put T3 to T
         T=zeros(12);
         for j=1:3:12
        T(j:j+2,j:j+2)=T(j:j+2,j:j+2)+T3;
         end
         
         %stiffness matrix for one element
         b=(elements_beam(i,5)*E)/(2*a);
         c=(3*E*Iz)/(2*a^3);
         d=(3*E*Iz)/(2*a^2);
         e=(3*E*Iy)/(2*a^3);
         f=(3*E*Iy)/(2*a^2);
         g=(G*J)/(2*a);
         h=(E*Iy)/a;
         k=(E*Iz)/a;
         
         k_element_t=[ b 0 0 0 0 0 -b 0 0 0 0 0 ; 0 c 0 0 0 d 0 -c 0 0 0 d ; 
        0 0 e 0 -f 0 0 0 -e 0 -f 0 ; 0 0 0 g 0 0 0 0 0 -g 0 0 ;0 0 -f 0 2*h 0 0 0 f 0 h 0; 
        0 d 0 0 0 2*k 0 -d 0 0 0 k ; -b 0 0 0 0 0 b 0 0 0 0 0;0 -c 0 0 0 -d 0 c 0 0 0 -d ; 
        0 0 -e 0 f 0 0 0 e 0 f 0; 0 0 0 -g 0 0 0 0 0 g 0 0 ; 0 0 -f 0 h 0 0 0 f 0 2*h 0;
        0 d 0 0 0 k 0 -d 0 0 0 2*k] ;
    
    ke=T'*k_element_t*T;
    
    %input the stifness matrix of its element in total stiffness K matrix
    %at the appropriate place according to nodes of element and direction
    %the series for its node to the appropriate place in K: j=1+(i-1)*6
    
    j=6*elements_beam(i,2)-5; %calculates the 1st cell of the first node of i element in GLOBAL Kstiff
    k=6*elements_beam(i,3)-5;  %calculates the 1st cell of the second node of i element in GLOBAL Kstiff
    K_stiff1(j:j+5,j:j+5)=K_stiff1(j:j+5,j:j+5)+ke(1:6,1:6); %places the 1st nodes x y ö èx èy èz
    K_stiff1(j:j+5,k:k+5)=K_stiff1(j:j+5,k:k+5)+ke(1:6,7:12); %each Ke : 12x12 matrix (6dof for 1st and 2nd node of each element)
    K_stiff1(k:k+5,j:j+5)=K_stiff1(k:k+5,j:j+5)+ke(7:12,1:6);
    K_stiff1(k:k+5,k:k+5)=K_stiff1(k:k+5,k:k+5)+ke(7:12,7:12); %places the 2nd nodes x y ö Èx èy èz
    
    else
        
        elements_beam_length(i,1)=norm(nodes_beam(elements_beam(i,3),:)-nodes_beam(elements_beam(i,2),:)); % distance between: (2nd node's x of element i)-(1st node's x of element i
        %direction x=l , sina, cosa etc
        l=(nodes_beam(elements_beam(i,3),2)-nodes_beam(elements_beam(i,2),2))/elements_beam_length(i,1);
        %direction y=m
        m=(nodes_beam(elements_beam(i,3),3)-nodes_beam(elements_beam(i,2),3))/elements_beam_length(i,1);
        %direction z=n
        n=(nodes_beam(elements_beam(i,3),4)-nodes_beam(elements_beam(i,2),4))/elements_beam_length(i,1);
        
        %stiffnes matrix for one element
         ke=(E*elements_beam(i,5)/elements_beam_length(i,1))*[l^2 l*m l*n 0 0 0 -l^2 -l*m -l*n 0 0 0 ;
        l*m m^2 m*n 0 0 0 -l*m -m^2 -m*n 0 0 0 ;
        l*n m*n n^2 0 0 0 -l*n -m*n -n^2 0 0 0 ;
        0 0 0 1 0 0 0 0 0 0 0 0; 
        0 0 0 0 1 0 0 0 0 0 0 0; 
        0 0 0 0 0 1 0 0 0 0 0 0;
       -l^2 -l*m -l*n 0 0 0 l^2 l*m l*n 0 0 0;
       -l*m -m^2 -m*n 0 0 0 l*m m^2 m*n 0 0 0;
       -l*n -m*n -n^2 0 0 0 l*n m*n n^2 0 0 0;
        0 0 0 0 0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0 0 0 1 0;
        0 0 0 0 0 0 0 0 0 0 0 1];
    
    %input the stifness matrix of its element in total stiffness K matrix
    %at the appropriate place according to nodes of element and direction
    %the series for its node to the appropriate place in K: j=1+(i-1)*6
    j=6*elements_beam(i,2)-5; %calculates the 1st cell of the first node of i element in GLOBAL Kstiff
    k=6*elements_beam(i,3)-5;  %calculates the 1st cell of the second node of i element in GLOBAL Kstiff
    K_stiff2(j:j+5,j:j+5)=K_stiff2(j:j+5,j:j+5)+ke(1:6,1:6); %places the 1st nodes x y ö èx èy èz
    K_stiff2(j:j+5,k:k+5)=K_stiff2(j:j+5,k:k+5)+ke(1:6,7:12); %each Ke : 12x12 matrix (6dof for 1st and 2nd node of each element)
    K_stiff2(k:k+5,j:j+5)=K_stiff2(k:k+5,j:j+5)+ke(7:12,1:6);
    K_stiff2(k:k+5,k:k+5)=K_stiff2(k:k+5,k:k+5)+ke(7:12,7:12); %places the 2nd nodes x y ö Èx èy èz
    end
end
    
%TOTAL STIFFNESS MATRIX
K_stiff=K_stiff1+K_stiff2;

%change K_stiff to K matrix according to Boundary Conditions
%Penalty Method
%where there is BC, multiply K(i,j) with large number
K=K_stiff; 
for i=1:size(bc_nodes_beam,1)
    if bc_nodes_beam(i,2)==1
        K(6*i-5,6*i-5)=K(6*i-5,6*i-5)+10^9;
    end
    if bc_nodes_beam(i,3)==1
        K(6*i-4,6*i-4)=K(6*i-1,6*i-4)+10^9;
    end
    if bc_nodes_beam(i,4)==1
        K(6*i-3,6*i-3)=K(6*i-3,6*i-3)+10^9;
    end
    if bc_nodes_beam(i,5)==1
         K(6*i-2,6*i-2)=K(6*i-2,6*i-2)+10^9;
    end
    if bc_nodes_beam(i,6)==1
         K(6*i-1,6*i-1)=K(6*i-1,6*i-1)+10^9;
    end
    if bc_nodes_beam(i,7)==1
         K(6*i,6*i)=K(6*i,6*i)+10^9;
    end
end

%create matrix with the forces on nodes with dimensions according to F=Ku
F=reshape(force_beam(1:32,2:7)',[],1);

%check if K is invertible
D=det(K);
if D==0
    msgbox('Error! Stiffness Matrix not invertible! Check Boundary Conditions')
    error('Error! K not invertible')
end

%Solve F=K*u
u=inv(K)*F; %nodes displacements
%Calculate Reaction Forces
Rs=K_stiff*u; 


%convert U to Matrix for 6dof
U=vec2mat(u,6);


%convert Reactions to Matrix for 6dof
Reactions=vec2mat(Rs,6);


%define strain matrix
strain=zeros(size(elements_beam,1),1);
stress=zeros(size(elements_beam,1),1);


%calculate new coordinations for nodes
for i=1:(size(nodes_beam,1)-1)
    nodes_newcord(i,1)=nodes_beam(i,1);
    nodes_newcord(i,2)=nodes_beam(i,2)+U(i,1);
    nodes_newcord(i,3)=nodes_beam(i,3)+U(i,2);
    nodes_newcord(i,4)=nodes_beam(i,4)+U(i,3);
    nodes_newcord(i,5)=nodes_beam(i,5)+U(i,4);
    nodes_newcord(i,6)=nodes_beam(i,6)+U(i,5);
    nodes_newcord(i,7)=nodes_beam(i,7)+U(i,6);
end

%calculate new length of each element
for i=1:size(elements_beam,1)
    elements_beam_newlength(i,1)=norm(nodes_newcord(elements_beam(i,3),2:4)-nodes_newcord(elements_beam(i,2),2:4));
    %calculate strains of its element
    strain(i,1)=(elements_beam_newlength(i,1)-elements_beam_length(i,1))/(elements_beam_length(i,1));
    stress(i,1)=210000*strain(i,1);
end


%calculate total reaction on each node
Rtotal(:,1)=sqrt(Reactions(:,1).^2+Reactions(:,2).^2+Reactions(:,3).^2+Reactions(:,4).^2+Reactions(:,5).^2+Reactions(:,6).^2);
Reactions=[Reactions,Rtotal];


%create txt file for nodes first coordinates
fID=fopen('AEM6101_Nodes_ger_beam.txt','w');
for i=1:size(nodes_beam,1)
        fprintf(fID,'%3d,   %6d,     %6d,    %6d,    %6d,      %6d,      %6d\r\n',nodes_beam(i,1),nodes_beam(i,2),nodes_beam(i,3),nodes_beam(i,4),nodes_beam(i,5),nodes_beam(i,6),nodes_beam(i,7));
end
fclose(fID);


%create txt file for nodes new coordinates
fID=fopen('AEM6101_NodesNewCord_ger_beam.txt','w');
for i=1:size(nodes_newcord,1)
    fprintf(fID,'%3d,    %15d,     %15d,     %15d,    %15d,      %15d,    %15d\r\n',nodes_newcord(i,1),nodes_newcord(i,2),nodes_newcord(i,3),nodes_newcord(i,4),nodes_newcord(i,5),nodes_newcord(i,6),nodes_newcord(i,7));
end
fclose(fID);


%create txt file for boundary conditions
fID=fopen('AEM6101_BC_ger_beam.txt','w');
for i=1:size(bc_nodes_beam,1)
    fprintf(fID,'%3d,       %5d,      %5d,       %5d,     %5d,     %5d,    %5d \r\n',bc_nodes_beam(i,1),bc_nodes_beam(i,2),bc_nodes_beam(i,3),bc_nodes_beam(i,4),bc_nodes_beam(i,5),bc_nodes_beam(i,6),bc_nodes_beam(i,7));
end
fclose(fID);




%create txt file for forces
fID=fopen('AEM6101_Forces_ger_beam.txt','w');
for i=1:size(force_beam,1)
    fprintf(fID,'%3d,   %5d,    %5d,     %5d,     %5d,    %5d,    %5d \r\n',force_beam(i,1),force_beam(i,2),force_beam(i,3),force_beam(i,4),force_beam(i,5),force_beam(i,6),force_beam(i,7));
end
fclose(fID);


%create txt file for elements
fID=fopen('AEM6101_Elements_ger_beam.txt','w');
for i=1:size(elements_beam,1)
    fprintf(fID,'%5d,   %4d,    %4d,   %8.4d,    %8d \r\n',elements_beam(i,1),elements_beam(i,2),elements_beam(i,3),elements_beam(i,4),elements_beam(i,5));
end
fclose(fID);

%create txt file for elements length
fID=fopen('AEM6101_Elementslength_ger_beam.txt','w');
for i=1:size(elements_beam_length,1)
    fprintf(fID,'%15d \r\n',elements_beam_length(i,1));
end
fclose(fID);

%create txt file for elements new length
fID=fopen('AEM6101_ElementsNewlength_ger_beam.txt','w');
for i=1:size(elements_beam_newlength,1)
    fprintf(fID,'%15d \r\n',elements_beam_newlength(i,1));
end
fclose(fID);

%create txt file for Strains 
fID=fopen('AEM6101_Strains_ger_beam.txt','w');
for i=1:size(strain,1)
    fprintf(fID,'%12.4d \r\n',strain(i,1));
end
fclose(fID);


%create txt file for Displacements
fID=fopen('AEM6101_Displacement_ger_beam.txt','w');
for i=1:size(U,1)
    fprintf(fID,'%11.4d,     %11.4d,      %11.4d,     %11.4d,     %11.4d,    %11.4d\r\n',U(i,1),U(i,2),U(i,3),U(i,4),U(i,5),U(i,6));
end
fclose(fID);


%create Reactions txt file
fID=fopen('AEM6101_Reactions_ger_beam.txt','w');
for i=1:size(Reactions,1)
    fprintf(fID,'%15.2d,   %15.2d,    %15.2d,  %15.2d,    %15.2d,  %15.2d,   %15.2d\r\n',Reactions(i,1),Reactions(i,2),Reactions(i,3),Reactions(i,4),Reactions(i,5),Reactions(i,6),Reactions(i,7));
end
fclose(fID);


%create Stresses txt file
fID=fopen('AEM6101_Stresses_ger_beam.txt','w');
for i=1:size(stress,1)
    fprintf(fID,'%12.4d \n',stress(i,1));
end
fclose(fID);






msgbox('Congratulations! Sucess! Good job! Proceed to Post processor')
