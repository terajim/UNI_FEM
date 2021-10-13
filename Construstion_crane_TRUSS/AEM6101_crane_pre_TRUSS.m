% pre_processor geranos
%Terzis Dimitrios AEM6101

a=0; % variable for while
while a~=1
    AEM=inputdlg('Enter your AEM (4digit number)','AEM Input');
    AEM=str2double(AEM);
    if isnumeric(AEM)==1 && numel(num2str(AEM))==4 %if input according to specifications then it's accepted
        a=1;
    else
        waitfor(msgbox('Wrong Input! Use 4digit Number for AEM')); %if wrong input while is continued
    end
end
[F,E,A,B,Angle,A0,L,Ad,Ay]=AEM6101_data_calc(AEM); %function created for data according to AEM
Storeys=5;

clear nodes


 %algorithmos vasis
 nodes=zeros(2,4);
for i=1:2                              
        nodes(i,1)=i;
        nodes(i,2)=A;
        nodes(i,3)=(i-1)*L;
        nodes(i,4)=0;
end

for i=1:Storeys+1
  for j=(4*i-1):(4*i+2)
        if (j-(4*(i-1)+2))==1
           nodes(j,1)=j;
           nodes(j,2)=A+(L/sind(60))*cosd(Angle-30)+(i-1)*L*cosd(Angle);
           nodes(j,3)=0;
           nodes(j,4)=(L/sind(60))*sind(Angle-30)+(i-1)*L*sind(Angle);
       elseif (j-(4*(i-1)+2))==2
            nodes(j,1)=j;
            nodes(j,2)=A+(L/sind(60))*cosd(Angle-30)+(i-1)*L*cosd(Angle);
            nodes(j,3)=L;
            nodes(j,4)=(L/sind(60))*sind(Angle-30)+(i-1)*L*sind(Angle);
        elseif (j-(4*(i-1)+2))==3
            nodes(j,1)=j;
            nodes(j,2)=A+(L/sind(60))*cosd(Angle+30)+(i-1)*L*cosd(Angle);
            nodes(j,3)=L;
            nodes(j,4)=(L/sind(60))*sind(Angle+30)+(i-1)*L*sind(Angle);
        else
            nodes(j,1)=j;
            nodes(j,2)=A+(L/sind(60))*cosd(Angle+30)+(i-1)*L*cosd(Angle);
            nodes(j,3)=0;
            nodes(j,4)=(L/sind(60))*sind(Angle+30)+(i-1)*L*sind(Angle);
            
        end
   end
end
 
 for i=1:3
     if i==3
        nodes(size(nodes,1)+1,1)=size(nodes,1)+1;
        nodes(size(nodes,1),2)=A+(L/sind(60))*cosd(Angle-30)+(Storeys+1)*L*cosd(Angle)+L;
        nodes(size(nodes,1),3)=L/2;
        nodes(size(nodes,1),4)=(L/sind(60))*sind(Angle-30)+(Storeys+1)*L*sind(Angle)-L/2;  
     elseif i==2
        nodes(size(nodes,1)+1,1)=size(nodes,1)+1;
        nodes(size(nodes,1),2)=A+(L/sind(60))*cosd(Angle-30)+(Storeys+1)*L*cosd(Angle);
        nodes(size(nodes,1),3)=L;
        nodes(size(nodes,1),4)=(L/sind(60))*sind(Angle-30)+(Storeys+1)*L*sind(Angle);
     else                  
        nodes(size(nodes,1)+1,1)=size(nodes,1)+1;
        nodes(size(nodes,1),2)=A+(L/sind(60))*cosd(Angle-30)+(Storeys+1)*L*cosd(Angle);
        nodes(size(nodes,1),3)=0;
        nodes(size(nodes,1),4)=(L/sind(60))*sind(Angle-30)+(Storeys+1)*L*sind(Angle);
     end
 end

 k=0;
 for i=1:3
        nodes(size(nodes,1)+1,1)=size(nodes,1)+1;
        nodes(size(nodes,1),2)=0;
        nodes(size(nodes,1),3)=0+k;
        nodes(size(nodes,1),4)=B;
        k=k+L/2;
 end


 
 
%Element Creation
clear elements
elements=ones(14+Storeys*20+4,3);
ID='TRUSS';
%elements vasis
elements(1,1)=1;
elements(1,2)=nodes(1,1);
elements(1,3)=nodes(2,1);
for i=2:5
    if i==2 || i==5
        elements(i,2)=nodes(1,1);
    else
        elements(i,2)=nodes(2,1);
    end
    elements(i,1)=i;
    elements(i,3)=nodes(i+1,1);
        
end

%perimetrika elements

for i=1:(Storeys+1)
    for j=(4*i+2):(4*i+5)
        if mod(j,4)==1
           elements(j,1)=j;
           elements(j,2)=nodes(j-3,1);
           elements(j,3)=nodes(j-6,1);
        else
            elements(j,1)=j;
            elements(j,2)=nodes(j-3,1);
            elements(j,3)=nodes(j-2,1);
        end
    end
end


%katheta elements (ektos tis korifis)
k=3;
for i=1:(Storeys)
    for j=(4*(Storeys+1)+4*i+2):(4*(Storeys+1)+4*i+5)
        elements(j,1)=j;
        elements(j,2)=nodes(k,1);
        elements(j,3)=nodes(k+4,1);
        k=k+1;
    end
end

%elements korifis


Logic1=true;
Logic2=true;
m=6;
v=3;
for j=(8*Storeys+10):8*Storeys+15
    elements(j,1)=j;
    elements(j,2)=nodes(size(nodes,1)-m-3,1);
    if Logic1==false
        elements(j,3)=nodes(size(nodes,1)-3);
    elseif Logic2==false
        elements(j,3)=nodes(size(nodes,1)-m-3+v);
        v=v-2;
    else
        elements(j,3)=nodes(size(nodes,1)-m-3+4);
    end
    m=m-1;
    if m==4
        if Logic1==false
            Logic2=false;
            Logic1=true;
            m=4;
        else
        Logic1=false;
        m=6;
        end
    end
end
n=5;
for j=(8*Storeys+16):8*Storeys+17
     elements(j,1)=j;
     elements(j,2)=nodes(size(nodes,1)-n,1);
     elements(j,3)=nodes(size(nodes,1)-3,1);
     n=n-1;
end
elements(8*Storeys+18,1)=8*Storeys+18;
elements(8*Storeys+18,2)=nodes(size(nodes,1)-5,1);
elements(8*Storeys+18,3)=nodes(size(nodes,1)-4,1);



%diagwnia elements
%vasis
o=4;
Logic=true;
for i=1:2
    for j=(8*Storeys+18+2*i-1):8*Storeys+18+2*i
        elements(j,1)=j;
        elements(j,2)=nodes(i,1);
        elements(j,3)=nodes(o,1);
        if Logic==false
            o=o+3;
        else
            o=o+1;
        end
    end
        o=3;
        Logic=false;
end

%ola ta upoloipa 

temp=3;
for i=1:Storeys
    orio1=temp+4;
    orio2=temp+6;
    for j=(8*Storeys+22+10*i-9):8*Storeys+22+10*i
        elements(j,1)=j;
        if temp<orio1
            elements(j,2)=nodes(temp,1);
            if mod(j,10)==3
                elements(j,3)=nodes(temp+7,1);
            else
                elements(j,3)=nodes(temp+3,1);
            end
        elseif temp<orio2
            elements(j,2)=nodes(temp-4,1);
            elements(j,3)=nodes(temp-2,1);
        else
           elements(j,2)=nodes(temp-2,1);
           if mod(j,10)==9
                elements(j,3)=nodes(temp-3,1); 
           else
                elements(j,3)=nodes(temp-7,1);
           end
        end
        temp=temp+1;
    end
    temp=temp-6;
end
             
%diagwnia korifis
j=Storeys*4+3;
l=j;
k=5;
Logic=false;
a=0;
for i=1:4
    while k>=4 && Logic==false
        elements(size(elements,1)-k,1)=size(elements,1)-k;
        elements(size(elements,1)-k,2)=nodes(j,1);
        elements(size(elements,1)-k,3)=nodes((j+k-1-a+1),1);
        k=k-1;
        j=j+1;
        a=a+1;
    end
    elements(size(elements,1)-k,1)=size(elements,1)-k;
    elements(size(elements,1)-k,2)=nodes(l,1);
    elements(size(elements,1)-k,3)=nodes(l+2,1); 
    l=l+1;
    Logic=true;
    k=k-1;
end

%sxoinia
elements(119,1)=119;
elements(119,2)=30;
elements(119,3)=14;
elements(120,1)=120;
elements(120,2)=31;
elements(120,3)=25;
elements(121,1)=121;
elements(121,2)=31;
elements(121,3)=26;
elements(122,1)=122;
elements(122,2)=32;
elements(122,3)=13;





%element properties
if strcmp(ID,'TRUSS')
    Properties=ones(size(elements,1),3);
    for i=1:8*Storeys+18
        Properties(i,1)=i;
        Properties(i,2)=Ay;
        Properties(i,3)=E;
    end
    for i=8*Storeys+19:size(elements,1)
        Properties(i,1)=i;
        Properties(i,2)=Ad;
        Properties(i,3)=E;
    end
else
    Properties=ones(4*Storeys+18,3);
    for i=1:8*Storeys+18
        Properties(i,1)=i;
        Properties(i,2)=Ay;
        Properties(i,3)=E;
    end
end

    
q1=1;
 bc_nodes=zeros(size(nodes,1),7); %boundary conditions matrix, 1 row for each node , 1-3 columns x y z , 4-6 columns èx èy èz 
 force=zeros(size(nodes,1),7);
 for i=1:size(nodes,1)
        bc_nodes(i,1)=i; %  Node ID
        force(i,1)=i;  %force Node ID
 end
 
 % ask user if wanna change boundary Contitions
 answer=listdlg('PromptString','Use Default boundary conditions?','SelectionMode','single','ListString',{'Yes','No'});
 
                %%%%-----------geometry pop up--------------%%%%
    fig = figure;
    
    str=ones(size(nodes,1),1);
    for j=1:size(nodes,1)
        str(j,1)=j;
        scatter3(nodes(j,2),nodes(j,3),nodes(j,4),'red','filled')
        text(nodes(j,2),nodes(j,3),nodes(j,4),{'node',num2str(str(j,1))})
        hold on
    end
    title('Geranos') %title of plots
    xlabel('x-axis [mm]')   %specify units of measurement
    ylabel('y-axis [mm]')   
    zlabel('z-axis [mm]')
    for i=1:8*Storeys+18
        line(nodes(elements(i,2:3),2),nodes(elements(i,2:3),3),nodes(elements(i,2:3),4),'LineWidth',3)
    end 
    for i=8*Storeys+19:size(elements,1)-3
        line(nodes(elements(i,2:3),2),nodes(elements(i,2:3),3),nodes(elements(i,2:3),4),'LineWidth',1)
    end 
    for i=size(elements,1)-3:size(elements,1)
        line(nodes(elements(i,2:3),2),nodes(elements(i,2:3),3),nodes(elements(i,2:3),4),'Color','b','LineWidth',1)
    end    
    axis equal
    camproj('orthographic')
    view([-8 -15 1])
    
                        %%%%-------------------------%%%%
 
 if answer==2
               
           %%%%-------------Boundary Conditions creation------------%%%%
    
    while q1==1
         list=listdlg('PromptString','Please select a Boundary Condition Style','SelectionMode','single','ListString',{'x,y,z','x,y','x,z','y,z','x','y','z'});
         hold on
         if list==1
             hold on
             shg
             dcm_obj = datacursormode(fig);
             set(dcm_obj,'DisplayStyle','window',...
             'SnapToDataVertex','on','Enable','on')
             camproj('orthographic')
             view([-8 -14 1])
             msgbox('Please select a node to lock its x,y,z')
             waitforbuttonpress
             c_info = getCursorInfo(dcm_obj);
             set(c_info.Target,'MarkerFaceColor','magenta')
             for i=1:size(nodes,1)
                if sum(c_info.Position==nodes(i,2:4))==3
                    bc_nodes(i,2:4)=1;
                end
             end
             datacursormode off
             
             q1=listdlg('PromptString','Do you wish to add more?','SelectionMode','single','ListString',{'Yes','No'});                
         elseif list==2
             hold on
             shg
             dcm_obj = datacursormode(fig);
             set(dcm_obj,'DisplayStyle','window',...
             'SnapToDataVertex','on','Enable','on')
             camproj('orthographic')
             view([-8 -14 1])
             msgbox('Please select a node to lock its x,y')
             waitforbuttonpress
             c_info = getCursorInfo(dcm_obj);
             set(c_info.Target,'MarkerFaceColor','yellow')
             for i=1:size(nodes,1)
                if sum(c_info.Position==nodes(i,2:4))==3
                    bc_nodes(i,2:3)=1;
                end
             end
             datacursormode off
             q1=listdlg('PromptString','Do you wish to add more?','SelectionMode','single','ListString',{'Yes','No'}); 
         elseif list==3
             hold on
             shg
             dcm_obj = datacursormode(fig);
             set(dcm_obj,'DisplayStyle','window',...
             'SnapToDataVertex','on','Enable','on')
             camproj('orthographic')
             view([-8 -14 1])
             msgbox('Please select a node to lock its x,z')
             waitforbuttonpress
             c_info = getCursorInfo(dcm_obj);
             set(c_info.Target,'MarkerFaceColor','yellow')
             for i=1:size(nodes,1)
                if sum(c_info.Position==nodes(i,2:4))==3
                    bc_nodes(i,2)=1;
                    bc_nodes(i,4)=1;
                end
             end
             datacursormode off
             q1=listdlg('PromptString','Do you wish to add more?','SelectionMode','single','ListString',{'Yes','No'}); 
         elseif list==4
             hold on
             shg
             dcm_obj = datacursormode(fig);
             set(dcm_obj,'DisplayStyle','window',...
             'SnapToDataVertex','on','Enable','on')
             camproj('orthographic')
             view([-8 -14 1])
             msgbox('Please select a node to lock its y,z')
             waitforbuttonpress
             c_info = getCursorInfo(dcm_obj);
             set(c_info.Target,'MarkerFaceColor','yellow')
             for i=1:size(nodes,1)
                if sum(c_info.Position==nodes(i,2:4))==3
                    bc_nodes(i,3:4)=1;
                end
             end
             datacursormode off
             q1=listdlg('PromptString','Do you wish to add more?','SelectionMode','single','ListString',{'Yes','No'}); 
         elseif list==5
             hold on
             shg
             dcm_obj = datacursormode(fig);
             set(dcm_obj,'DisplayStyle','window',...
             'SnapToDataVertex','on','Enable','on')
             camproj('orthographic')
             view([-8 -14 1])
             msgbox('Please select a node to lock its x')
             waitforbuttonpress
             c_info = getCursorInfo(dcm_obj);
             set(c_info.Target,'MarkerFaceColor','black')
             for i=1:size(nodes,1)
                if sum(c_info.Position==nodes(i,2:4))==3
                    bc_nodes(i,2)=1;
                end
             end
             datacursormode off
             q1=listdlg('PromptString','Do you wish to add more?','SelectionMode','single','ListString',{'Yes','No'}); 
         elseif list==6
             hold on
             shg
             dcm_obj = datacursormode(fig);
             set(dcm_obj,'DisplayStyle','window',...
             'SnapToDataVertex','on','Enable','on')
             msgbox('Please select a node to lock its y')
             camproj('orthographic')
             view([-8 -14 1])
             waitforbuttonpress
             c_info = getCursorInfo(dcm_obj);
             set(c_info.Target,'MarkerFaceColor','black')
             for i=1:size(nodes,1)
                if sum(c_info.Position==nodes(i,2:4))==3
                    bc_nodes(i,3)=1;
                end
             end
             datacursormode off
             q1=listdlg('PromptString','Do you wish to add more?','SelectionMode','single','ListString',{'Yes','No'}); 
         else
             hold on
             shg
             dcm_obj = datacursormode(fig);
             set(dcm_obj,'DisplayStyle','window',...
             'SnapToDataVertex','on','Enable','on')
             camproj('orthographic')
             view([-8 -14 1])
             msgbox('Please select a node to lock its z')
             waitforbuttonpress
             c_info = getCursorInfo(dcm_obj);
             set(c_info.Target,'MarkerFaceColor','black')
             for i=1:size(nodes,1)
                if sum(c_info.Position==nodes(i,2:4))==3
                    bc_nodes(i,4)=1;
                end
             end
             datacursormode off
             q1=listdlg('PromptString','Do you wish to add more?','SelectionMode','single','ListString',{'Yes','No'}); 
         end
         
    end
    
 else
    bc_nodes=zeros(size(nodes,1),7); %boundary conditions matrix, 1 row for each node , 1-3 columns x y z , 4-6 columns èx èy èz 
    for i=1:size(nodes,1)
        bc_nodes(i,1)=i;
    end
    hold on
    bc_nodes(1,2:4)=1;
    scatter3(nodes(1,2),nodes(1,3),nodes(1,4),'magenta','filled')
    bc_nodes(2,2)=1;
    bc_nodes(2,4)=1;
    scatter3(nodes(2,2),nodes(2,3),nodes(2,4),'yellow','filled')
    bc_nodes(size(nodes,1)-2:size(nodes,1),2:4)=1;
    scatter3(nodes(size(nodes,1)-2:size(nodes,1),2),nodes(size(nodes,1)-2:size(nodes,1),3),nodes(size(nodes,1)-2:size(nodes,1),4),'magenta','filled')
 end
           
    
    hold off
    
   
    %%%%-------------Force position and value------------%%%%
                     
    
    forceq=listdlg('PromptString','Change the force node?','SelectionMode','single','ListString',{'Yes','No'}); 
    if forceq==1
       k=3;
       while k~=1 
           hold on
           shg
           dcm_obj = datacursormode(fig);
           set(dcm_obj,'DisplayStyle','window',...
           'SnapToDataVertex','on','Enable','on')
           msgbox('Please set the force node')
           waitforbuttonpress
           c_info = getCursorInfo(dcm_obj);
           set(c_info.Target,'MarkerFaceColor','green')
           force_ans=inputdlg({'Force (N)= ','Coordinate: 2 for x, 3 for y, 4 for z , 5 for èx,6 for èy,7 for èz'},'Input',[1 35],{'-20000','4'});
           force_ans=str2double(force_ans);
           if isnumeric(force_ans)==1 && force_ans(2)>=2 && force_ans(2)<=7
               k=1;
           else
               waitfor(msgbox('Wrong Input!'));
           end
       end
       for i=1:size(nodes,1)
           if sum(c_info.Position==nodes(i,2:4))==3
               force(i,force_ans(2))=force_ans(1);
           end
       end
    else
         hold on
         force(29,4)=-20000;
         scatter3(nodes(29,2),nodes(29,3),nodes(29,4),'green','filled')
         hold off
    end
    hold on
    h1=scatter(NaN,NaN,'magenta','LineWidth',3);
    h2=scatter(NaN,NaN,'yellow','LineWidth',3);
    h3=scatter(NaN,NaN,'black','LineWidth',3);
    h4=scatter(NaN,NaN,'green','LineWidth',3);
    h5=scatter(NaN,NaN,'red','LineWidth',3);
    h6=scatter(NaN,NaN,'blue','LineWidth',3);
    legend([h1,h2,h3,h4,h5,h6],'0 dof','1 dof','2 dof ','Force','Nodes','Elements')
    hold off
   
    
    
    
    %create txt file for nodes
    fID=fopen('AEM6101_nodes_ger_truss.txt','w');
    for i=1:size(nodes,1)
        fprintf(fID,'%d,    %4d,     %4d,     %4d\r\n',nodes(i,1),nodes(i,2),nodes(i,3),nodes(i,4));
    end
    fclose(fID);
    
    %create txt file for force applied
    fID=fopen('AEM6101_force_ger_truss.txt','w');
    for i=1:size(force,1)
        fprintf(fID,'%d,    %4d,      %4d,      %4d\r\n',force(i,1),force(i,2),force(i,3),force(i,4));
    end
    fclose(fID);
    
    %create txt file for boundary conditions
    fID=fopen('AEM6101_BC_ger_truss.txt','w');
    for i=1:size(bc_nodes,1)
        fprintf(fID,'%d,    %d,      %d,      %d\r\n',bc_nodes(i,1),bc_nodes(i,2),bc_nodes(i,3),bc_nodes(i,4));
    end
    fclose(fID);
    
    %create txt file for elements
    fID=fopen('AEM6101_elements_ger_truss.txt','w');
    for i=1:size(elements,1)
        fprintf(fID,'%d,    %d,     %d,     %d,     %d\r\n',elements(i,1), elements(i,2), elements(i,3), Properties(i,2), Properties(i,3));
    end
    fclose(fID);


 %create txt file for nodes
    fID=fopen('AEM6101_nodes_ger_truss.txt','w');
    for i=1:size(nodes,1)
        fprintf(fID,'%d,    %4d,     %4d,     %4d\r\n',nodes(i,1),nodes(i,2),nodes(i,3),nodes(i,4));
    end
    fclose(fID);