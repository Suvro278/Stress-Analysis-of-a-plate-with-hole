% Please note that this file is written specificallyf or cantilever beam
% problem

% bc_node_x = this vector stores the nodes where x dof is fixed
% bc_node_y = this vector stores the nodes where y dof is fixed


bc_node_x = [] ;
bc_node_y = [] ;

number_of_nodes_with_x_bc =  0 ;
number_of_nodes_with_y_bc =  0 ;

% Get the nodes which are lying on the Y axis as this end of the cantilever
% beam is fixed.
for i = 1:nno
    if abs(Xn(i,1)) <= eps % if x coordinate is zero
        bc_node_x = [bc_node_x ; i];
        %bc_node_y=[bc_node_y ; 1];
        number_of_nodes_with_x_bc = number_of_nodes_with_x_bc + 1 ;
        %number_of_nodes_with_y_bc = number_of_nodes_with_y_bc + 1 ;
    end
    if  abs(Xn(i,2)) <= 0.0001
        bc_node_y=[bc_node_y ; i];
        %bc_node_x = [bc_node_x ; ];
        number_of_nodes_with_y_bc = number_of_nodes_with_y_bc + 1 ;
    end
end


% Form the boundary condition array. ir contains the nodes where the
% displacements are unknown. To get this array we first have to get the 
% dofs in the x and y directions. These are stored in arrays ir_x and ir_y.
% For your own problem you have to 

% ir = zeros(number_of_nodes_with_x_bc+number_of_nodes_with_y_bc, 1) ;

% Form the boundary condition array in x direction

ir_x =  zeros(number_of_nodes_with_x_bc,1) ;

for i = 1: number_of_nodes_with_x_bc
    ir_x(i,1) = 2 * bc_node_x(i,1)-1 ; 
end

% Form the boundary condition array in y direction
ir_y =  zeros(number_of_nodes_with_y_bc,1) ;
for i = 1: number_of_nodes_with_y_bc
    ir_y(i,1) = 2 * bc_node_y(i,1); 
end


% Final boundary condition array. This will be used to apply the boundary
% condition on the stiffness matrix and external force vector

% This array constains all the dofs which are fixed i.e. have no displacement 

ir_disp = [ir_x ; ir_y];

% Now using setdiff function of MATLAB we get the list of dofs where the
% displacement is unknown
ir = setdiff(1:ndof,ir_disp)' ;