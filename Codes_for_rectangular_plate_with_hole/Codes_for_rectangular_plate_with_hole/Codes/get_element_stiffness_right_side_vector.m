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

function [kel, rel, ie, jterm] = get_element_stiffness_right_side_vector(...
    con,ndoel,oI,xI,U,ngpv,xigv,I2,D) 
    
%
% Input
%
% con - connectivity for ith element
% nodel - number of dof per eleemnt
% oI - undefomred coordinate
% xI - deformed coordinate
% U - displacement 
% ngpv - number of Gauss point required in the bulk integration
% xigv - this array contains the Gauss point 

%
% Output
% 
% kel - element stiffness matrix
% rel - elemental external force vector
 % ie - degree of freedom vector
% jterm - error flag 


jterm = 0 ;

% Global Dofs
ie = [con(1)-1 con(1) con(2)-1 con(2) con(3)-1 con(3) con(4)-1 con(4)] ;

% Displacements
uI = U(ie) ; 

% Initialize 
rel  = zeros(ndoel,1) ;
kmat = zeros(ndoel,ndoel) ;
kgeo = zeros(ndoel,ndoel) ;
dN   = zeros(4,2) ;

kel  = zeros(ndoel,ndoel) ;

% Gauss Point Loop (Volume Integrals)
for gp = 1:ngpv

    % Gauss point coordinates, weight
    xi = xigv(gp,1) ; eta = xigv(gp,2) ; wg = xigv(gp,3) ;

    % Shape Functions
    N1 = ( 1 - xi )*( 1 - eta )/4 ;
    N2 = ( 1 + xi )*( 1 - eta )/4 ;
    N3 = ( 1 + xi )*( 1 + eta )/4 ;
    N4 = ( 1 - xi )*( 1 + eta )/4 ;
    N  = [N1*I2 N2*I2 N3*I2 N4*I2] ;
    
   % Shape Fct (xi,eta)-Derivatives
    dpN = [ -( 1-eta )  -( 1-xi ) ;
             ( 1-eta )  -( 1+xi ) ;
             ( 1+eta )   ( 1+xi ) ;
            -( 1+eta )   ( 1-xi ) ]/4 ;
    
    % Jacobian between Current Configuration and Master Element
    Jac = zeros(2,2) ;
    for j = 1:4
       Jac = Jac + xI(j,:)' * dpN(j,:) ;
    end
    detJ = det(Jac) ; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Jacobian between Reference Configuration and Master Element
    JacRef = zeros(2,2) ;
    for j = 1:4
       JacRef = JacRef + oI(j,:)' * dpN(j,:) ;
    end
    detJRef = det(JacRef) ;
    invJRef = inv(JacRef);       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Deformation Warning
    if detJ < 0       
        jterm = 1 ; % in case we have negative Jacobian then error is thrown
    end
    
    invJ = inv(Jac) ;
    
    % Shape Fct (x,y)- Spatial Derivatives
    for j = 1:4
       dN(j,:)  =  dpN(j,:) * invJ ;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Added on 27.07.2010 %%%%%%%%%%%%%%%%%%%%%%%
    % Shape Fct (x,y)- Spatial Derivatives
    for j = 1:4
       dNRef(j,:)  =  dpN(j,:) * invJRef ;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%% Added on 27.07.2010 %%%%%%%%%%%%%%%%%%%%%%%
    
    % B-Matrix           
    B = [ dN(1,1)        0  dN(2,1)        0  dN(3,1)        0  dN(4,1)        0 ;
                0  dN(1,2)        0  dN(2,2)        0  dN(3,2)        0  dN(4,2) ;
          dN(1,2)  dN(1,1)  dN(2,2)  dN(2,1)  dN(3,2)  dN(3,1)  dN(4,2)  dN(4,1) ] ;
      
    
    % Deformation Gradient
    invF = zeros(2,2) ;
    for j = 1:4
       invF = invF + oI(j,:)' * dN(j,:) ; 
    end
    Fgr  = inv(invF) ; % This is the formation gradient
    detF = det(Fgr) ;
    
    if detF < 0 
        jterm = 1 ;              
    end  
    
 
    % Material Stiffness
    kmat = kmat + B' * D * B * detJ * wg ;

    % Note that if you have distributed loads on some face then the
    % coresponding external force vector has to be considered here.

end % End of Gauss point loop

% Local Element Stiffness
kel = kmat ;

end