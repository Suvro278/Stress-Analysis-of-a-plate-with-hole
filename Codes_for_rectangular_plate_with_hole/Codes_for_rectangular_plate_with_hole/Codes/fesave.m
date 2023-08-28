
% --------> 1
% File for printing the initial coordinates
fileInitialCoordinate = 'Output/initial_coordinate.txt' ;
 % Opening the file for the first time
fidInitialCoordinate = fopen(fileInitialCoordinate,'w') ;
for i = 1:nno % run the loop over all the nodes
    fprintf(fidInitialCoordinate,'%20.15f \t %20.15f \n',Xn(i,1),Xn(i,2));
end
fclose(fidInitialCoordinate) ; % Close the file

% --------> 2
% File for printing the deformed  coordinates
fileDeformedCoordinate = 'Output/deformed_coordinate.txt' ;
% Opening the file for the first time
fidDeformedCoordinate = fopen(fileDeformedCoordinate,'w') ; 
for i = 1:nno
    fprintf(fidDeformedCoordinate,'%20.15f \t %20.15f \n',xn(i,1),xn(i,2));
end
fclose(fidDeformedCoordinate) ; % Close the file


% --------> 3
% File for printing the deformed  displacement
fileDeformedDisplacement = 'Output/deformed_displacement.txt' ;
% Opening the file for the first time
fidDeformedDisplacement = fopen(fileDeformedDisplacement,'w') ; 
for i = 1:nno % run the loop over all the nodes
    fprintf(fidDeformedDisplacement,'%20.15f \t %20.15f \n',Ux(i,1),Ux(i,2));
end
fclose(fidDeformedDisplacement) ; % Close the file


% --------> 4
% File for printing the displacement of the tip of the cantilever located 
% at the coordinates X = L Y = H
fileTipNodeCoordinate = 'Output/output_tip_node_coordinate.txt' ;
fidTipNodeCoordinate = fopen(fileTipNodeCoordinate,'w') ; % Opening the file for the first time
for i = 1:nno % run the loop over all the nodes
    fprintf(fidTipNodeCoordinate,'\n Initial coordinate of point with tip load is : % 20.15f,\t %20.15f  \n',Xn(end,1),Xn(end,2));
    fprintf(fidTipNodeCoordinate,'\n Deformed coordinate of point with tip load is : % 20.15f,\t %20.15f  \n',xn(end,1),xn(end,2));
    fprintf(fidTipNodeCoordinate,'\n Displacement of the point with tip load is : % 20.15f,\t %20.15f  \n',Ux(end,1),Ux(end,2));
end
fclose(fidTipNodeCoordinate) ;

% --------> 5
% Print the coordinates (undeformed and deformed) as well as displacement
% of the node at the end of the cantilever. Witht he way nodes are number
% in file preprocessor
%fprintf('\n Initial coordinate of point with tip load is : % 20.15f,\t %20.15f  \n',Xn(end,1),Xn(end,2));
%fprintf('\n Deformed coordinate of point with tip load is : % 20.15f,\t %20.15f  \n',xn(end,1),xn(end,2));

%fprintf('\n Displacement of the point with tip load is : % 20.15f,\t %20.15f  \n',Ux(end,1),Ux(end,2));