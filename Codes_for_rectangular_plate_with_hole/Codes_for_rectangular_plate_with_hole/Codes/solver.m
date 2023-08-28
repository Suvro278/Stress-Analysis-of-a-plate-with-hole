%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 INDIAN INSTITUTE OF TECHNOLOGY GUWAHATI                 %
%                  DEPARTMENT OF MECHANICAL ENGINEERING                   %
%                                                                         %
%                          2022-23 2ND SEMESTER                           %
%                                                                         %
%               ME 682 - NONLINEAR FINITE ELEMENT METHODS                 %
%                                                                         %
%                                                                         %
% Code initially developed by: Sachin Singh Gautam                        %
%                                                                         %
%                                                                         %
% Project 1: Due date 31.03.2023, Friday, 5 PM                            %
%                                                                         %
% The code is written for solving a finding the displacement, strains and % 
% stresses for a cantilever beam subjected to point load as shown below   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
I2 = eye(2) ; 
ksp = zeros(nksz,1) ; ndoel = ndoelo ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%                    Bulk Element Loop                                   %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kel_dummy = zeros(ndoel,ndoel,nel) ; rel_dummy = zeros(ndoel,nel) ;
jtermu = zeros(nel,1) ;
parfor i = 1:nel   % Run the loop over elements

   kel = zeros(ndoel,ndoel) ; rel = zeros(ndoel,1) ; 

   [kel, rel, ie, jtermue] = get_element_stiffness_right_side_vector(...
       2*CON(i,:),ndoel,Xn(CON(i,:),:),xn(CON(i,:),:),U,ngpv,xigv,I2,D) ; 

    kel_dummy(:,:,i) = kel ; 
    rel_dummy(:,i) = rel ; % This is external laod vector due to distributed load.
    jtermu(i,1) = jtermue ; % Store error flag
end

for i = 1:nel   
   con = 2*CON(i,:) ;
   ie = [con(1)-1 con(1) con(2)-1 con(2) con(3)-1 con(3) con(4)-1 con(4)] ;
   F_ext(ie,1) = F_ext(ie,1) + rel_dummy(:,i) ;
end

for i = 1:nel
  for j = 1:ndoel
   ksp(cspa_global(j,i):csp_global(j,i))  = kel_dummy(:,j,i) ;
  end
end

clear kel_dummy ;
clear rel_dummy ;

% Check if there was any error in the element loop then kindly throw the
% error message
if jterm == 0
    if sum(jtermu(:,1)) > 0
        jterm = 1 ;
        fprintf('\n\n Error in computation in element loop ..... exiting the newton iteration loop ...') ;
%             fprintf(fipf,'\n\n Error in computation in element loop ..... exiting the newton iteration loop ...') ; 
%             break ;        
    end
end

% Sparse Stiffness Matrix
K = sparse(isp,jsp,ksp) ;

clear ksp ;

FACE_CON=load("Input/FACE_CON.txt");
filename="Input\elemental_load_face_1.txt";
fid=fopen("Input\elemental_load_face_1.txt","w");
for i=(1+(nr-1)*na):[(nr*na)-0.5*na]
    fprintf(fid,"%d\t 3\n",i);
end
filename="Input\elemental_load_face_2.txt";
fid=fopen("Input\elemental_load_face_2.txt","w");
for i=(1+(nr-1)*na+0.5*na):[(nr*na)]
    fprintf(fid,"%d\t 3\n",i);
end

elemental_load_face_1=load("Input\elemental_load_face_1.txt");
elemental_load_face_2=load("Input\elemental_load_face_2.txt");
ta=10000000;
tb=20000000;
t1 = [0 ;1000000];
t2 = [1000000 ;0];
for i=1:(0.5*na)
    element_no=elemental_load_face_1(i,1);
    element_face=elemental_load_face_1(i,2);
    local_Con_face=FACE_CON(element_face,:);
    global_Con_face=CON(element_no,[local_Con_face]);
    ie_global_1=2*global_Con_face;
    ie_global=[ie_global_1(1)-1 ie_global_1(1) ie_global_1(2)-1 ie_global_1(2)];
    force_vector=zeros(4,1);
    
    

    for j=1:ngps
        a=Xn(i+1,1)-Xn(i,1);
        b=Xn(i+1,1)+Xn(i,1);
        %p=ta+((tb-ta)/Lx)*0.5*(b+(a*xigs(j,1)));
        %t1 = [0 ;p];
        N1=(1-xigs(j,1))/2;
        N2=(1+xigs(j,1))/2;
        N  = [N1 0 N2 0;0 N1 0 N2]; % size = [2x1]
        J=0.5*sqrt(((Xn(global_Con_face(1),1)-Xn(global_Con_face(2),1))^2)+((Xn(global_Con_face(1),2)-Xn(global_Con_face(2),2))^2));
        force_vector=force_vector+N'*t1*J*xigs(j,2);
    end  
     F_ext(ie_global)=F_ext(ie_global)+force_vector;
end
for i=1:0.5*na
    element_no=elemental_load_face_2(i,1);
    element_face=elemental_load_face_2(i,2);
    local_Con_face=FACE_CON(element_face,:);
    global_Con_face=CON(element_no,[local_Con_face]);
    ie_global_1=2*global_Con_face;
    ie_global=[ie_global_1(1)-1 ie_global_1(1) ie_global_1(2)-1 ie_global_1(2)];
    force_vector=zeros(4,1);
    

    for j=1:ngps
        a=Xn(i+1,2)-Xn(i,2);
        b=Xn(i+1,2)+Xn(i,2);
        %q=ta+((tb-ta)/Ly)*0.5*(b+(a*xigs(j,1)));
        %t2=[q;0];       
        N1=(1-xigs(j,1))/2;
        N2=(1+xigs(j,1))/2;
        N  = [N1 0 N2 0;0 N1 0 N2]; % size = [2x1]
        J=0.5*sqrt(((Xn(global_Con_face(1),1)-Xn(global_Con_face(2),1))^2)+((Xn(global_Con_face(1),2)-Xn(global_Con_face(2),2))^2));
        force_vector=force_vector+N'*t2*J*xigs(j,2);
    end  
     F_ext(ie_global)=F_ext(ie_global)+force_vector;
end



% Eliminate Constraints in K and R due to BC and Condensation
Kr = K(ir,ir) ;
fr = F_ext(ir,1); 

clear K ; clear f ;
   
% Get the sparse format K
Kr     = sparse(Kr) ;
  
% Solve for incremental displacement
DUr    =  Kr \ fr ;
   
clear Kr ;
   
DU     = zeros(ndof,1) ; % DU is the iterative displacement
DU(ir) = DUr ;  
U      = U + DU ; % Add incremental displacement to  previous displacement 
                  %to get the current total displacement
   
 % Current configuration
Ux = [ U(i1) U(i2) ];  % Rearragne the displacement in the same form as coordinate array
xn = Xn + Ux(:,1:2) ;