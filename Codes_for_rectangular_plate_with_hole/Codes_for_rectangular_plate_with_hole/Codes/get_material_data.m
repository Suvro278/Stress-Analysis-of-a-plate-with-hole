
% Bulk Material
E  = Eo/Lo^3    ; 

% Lames constant
mu = E/2/(1+nu) ; lam = 2*mu*nu/(1-2*nu) ; kappa = lam + 2*mu/3 ;

% Define the density  of the body (non-dimesionalzed)
rho = 0.5e-06;


if option_type_2D == 1 % plane strain
    
%     % Linear Material Tangent for plane strain case
%     Ip4 = [ones(2,2) zeros(2,1) ; zeros(1,2) 0 ] ;
%     Is4 = [ 2*eye(2) zeros(2,1) ; zeros(1,2) 1 ] ;
%     D   = lam*Ip4 + mu*Is4 ;

   % Linear Material Tangent for plane stress case
    term_1 = E/(((1.0 + nu))*(1.0 - 2 * nu * nu)) ;
    term_2 = (1 - nu) * term_1 ;
    term_3 = nu * term_1 ;
    term_4 = ((1 - 2*nu)/2) * term_1 ;
    
    D   = [ term_2   term_3          0
            term_3   term_2          0
                 0        0     term_4 ] ;             

end

if option_type_2D == 2 % plane stress

    % Linear Material Tangent for plane stress case
    term_1 = E/(1.0 - 2 * nu * nu) ;
    term_2 = nu * term_1 ;
    term_3 = (1.0 - nu)/2.d0 ;
    term_4 = term_1 * term_3 ;
    
    D   = [ term_1   term_2          0
            term_2   term_1          0
                 0        0     term_4 ] ;

end