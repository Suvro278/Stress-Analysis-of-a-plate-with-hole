% FEQUAD - Gaussian Quadature Points and Weights
% Sachin Singh Gautam - 05.03.23


% 1D Quadrature
if jgaus == 1
   gfa = 0.0;
   xig1(1,1) = gfa ;
   xig1(1,2) = 2.d0 ;
end

% 2 GP Quadrature
if jgaus == 2
  gfa = 1.d0/sqrt(3.d0) ;
  xig1(1,1) = -gfa ;
  xig1(1,2) = 1.d0 ;
  xig1(2,1) = gfa  ;
  xig1(2,2) = 1.d0 ; 
    
% 3 GP Quadrature
elseif jgaus == 3
  gfa = sqrt(3.d0/5.d0) ;
  xig1(1,1) = -gfa      ;
  xig1(1,2) = 5.d0/9.d0 ;
  xig1(2,1) = 0.d0      ;
  xig1(2,2) = 8.d0/9.d0 ;
  xig1(3,1) = gfa       ;
  xig1(3,2) = 5.d0/9.d0 ;   

% 5 GP Quadrature
elseif jgaus == 5
  xig1(1,1) = -0.906179845938664 ;
  xig1(1,2) =  0.236926885056189 ;
  xig1(2,1) = -0.538469310105683 ;
  xig1(2,2) =  0.478628670499366 ;
  xig1(3,1) =  0                 ;   
  xig1(3,2) =  0.568888888888889 ;  
  xig1(4,1) =  0.538469310105683 ;
  xig1(4,2) =  0.478628670499366 ;       
  xig1(5,1) =  0.906179845938664 ;
  xig1(5,2) =  0.236926885056189 ; 
  
% 10 GP Quadrature
elseif jgaus == 10
  xig1(1,1)  = -0.973906528517172 ;
  xig1(1,2)  =  0.066671344308688 ;
  xig1(2,1)  = -0.865063366688985 ;
  xig1(2,2)  =  0.149451349150581 ;
  xig1(3,1)  = -0.679409568299024 ;
  xig1(3,2)  =  0.219086362515982 ;   
  xig1(4,1)  = -0.433395394129247 ;
  xig1(4,2)  =  0.269266719309996 ;       
  xig1(5,1)  = -0.148874338981631 ;
  xig1(5,2)  =  0.295524224714753 ;
  xig1(6,1)  =  0.148874338981631 ;
  xig1(6,2)  =  0.295524224714753 ;
  xig1(7,1)  =  0.433395394129247 ;
  xig1(7,2)  =  0.269266719309996 ; 
  xig1(8,1)  =  0.679409568299024 ;
  xig1(8,2)  =  0.219086362515982 ;
  xig1(9,1)  =  0.865063366688985 ;
  xig1(9,2)  =  0.149451349150581 ;
  xig1(10,1) =  0.973906528517172 ;
  xig1(10,2) =  0.066671344308688 ;

else
%   %break  
%   for i = 1:jgaus
%     xig1(i,1) = i*2/jgaus-1/jgaus-1 ;
%     xig1(i,2) = 2/jgaus ;
%   end

[x,w]=lgwt(jgaus,-1,1) ;
xig1 = [x w ] ;

end
        

% 2D Quadrature
gp = 0 ;
for j = 1:jgaus
  for i = 1:jgaus
    gp = gp + 1 ;
    xig2(gp,1) = xig1(i,1) ;
    xig2(gp,2) = xig1(j,1) ;     
    xig2(gp,3) = xig1(i,2)*xig1(j,2) ;
  end
end