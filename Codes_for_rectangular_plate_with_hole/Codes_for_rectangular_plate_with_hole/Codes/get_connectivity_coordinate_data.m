

% Number of elements in y direction
nr=30;     %no of elements in axial direction
% Number of elements in x direction % <--- change this for your own problem
na=40;
% Number of Elements
nel=nr*na;
% Number of Nodes
nno = (nr+1)*(na+1) ;

% Number of dofs per element 
ndoel = 8 ;
 
% Number of dofs in the FE mesh, Number of dofs per element
ndof = 2*nno ; ndoelo = 8 ;

filename = 'Input/fe_data.txt' ;
fid = fopen(filename,'w') ;
fprintf(fid,'%g \t %g \t %g \t %g \t %g',nr,na,nel,ndoel,ndof);
fclose(fid);

% Element Connectivity
%
% Note that the below code works only for rectangular bodies. You can
% either write your own code or take input from comemrcial FE packages like
% ANSYS. In that case change the below code accordingly.

%CON = zeros(nel,4) ; % since 4 nodes per element
%iel = 0 ;
%for j = 1:ny
    %for i = 1:nx 
        %iel = iel + 1 ; 
        %CON(iel,:) = [ (j-1)*(nx+1)+i  (j-1)*(nx+1)+i+1  j*(nx+1)+i+1  j*(nx+1)+i ] ;    
    %end
%end
CON = zeros(nel,4) ; % since 4 nodes per element
iel = 0 ;
for j = 1:nr
    for i = 1:na 
        iel = iel + 1 ; 
        CON(iel,:) = [ (j-1)*(na+1)+i  (j-1)*(na+1)+i+1  j*(na+1)+i+1  j*(na+1)+i ] ;    
    end
end

filename = 'Input/connectivity.txt' ;
fid = fopen(filename,'w') ;
for iel = 1:nel
    if iel < nel
        fprintf(fid,'%g \t %g \t %g \t %g \n',CON(iel,:));   
    else
        fprintf(fid,'%g \t %g \t %g \t %g',CON(iel,:));
    end
end
fclose(fid);

X=zeros((na+1),2);
ino=0;
dtheta=0;
for i=1:(0.5*na+1)
    ino=ino+1;
    m=tan((pi/2)-dtheta);
    X(ino,:)=[Ly/(m) Ly];
    dtheta=dtheta+((pi/2)/na);
    
end
ino=(0.5*na+1);
theta=atan(Ly/Lx)-((pi/2)/na);
for i=1:(0.5*na)
    ino=ino+1;
    m=tan(theta);
    X(ino,:)=[Lx m*Lx];
    theta=theta-((pi/2)/na);
end

R=20;
dr=zeros((nr+1),1);
for i=1:size(X,1)
    d=sqrt(((X(i,1))^2)+((X(i,2))^2)); 
    dr(i)=(d-R)/nr;
end

dr;

Xn = zeros(nno,2) ;
ino=0;

for i=1:nr+1
    theta=0;
    dtheta=((pi/2)/na);
    for j=1:na+1
        ino=ino+1;
        Xn(ino,:)=[(R+(i-1)*dr(j))*sin(theta) (R+(i-1)*dr(j))*cos(theta)];
        theta=theta+dtheta;
    end
end




% Initial Node Coordinates
%Xn = zeros(nno,2) ; % 2 dof per node i.e. x coordinate and y coordinate.
%ino = 0 ;
%for j = 0:ny
    %for i = 0:nx 
        %ino = ino + 1 ;
        %Xn(ino,:) =  [ Xo + i*Lx/nx  Yo + j*Ly/ny ] ;
    %end
%end

filename = 'Input/coordinate.txt' ;
fid = fopen(filename,'w') ;
for ino = 1:nno
    if ino < nno
        fprintf(fid,'%20.13f \t %20.13f \n',Xn(ino,:));   
    else
        fprintf(fid,'%20.13f \t %20.13f',Xn(ino,:));   
    end
end
fclose(fid);


% Current node coordinates initialized as initial coordinates to start
% with.
xn = Xn ; % Xn contains the coordinates of the reference configuration all the time
          % xn contains the coordinates of the current configuration all
          % the time
