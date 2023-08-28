% FELP Principal Stresses 
% Roger Sauer - 14.10.09


% Principal Stress
R(:,7)  = (R(:,1)+R(:,2))/2 + sqrt( ((R(:,1)-R(:,2)).^2)/4 + R(:,4).^2 ) ;
R(:,8)  = (R(:,1)+R(:,2))/2 - sqrt( ((R(:,1)-R(:,2)).^2)/4 + R(:,4).^2 ) ;
    
% Invariants
R(:,9)  = R(:,1) + R(:,2) + R(:,3) ;
R(:,10) = R(:,7).*R(:,8) + R(:,8).*R(:,3) + R(:,3).*R(:,7) ;
R(:,11) = R(:,7).*R(:,8).*R(:,3) ;
    
% Von Mises
R(:,12) = sqrt(( (R(:,7)-R(:,8)).^2 + (R(:,8)-R(:,3)).^2 + (R(:,3)-R(:,7)).^2 )/2) ;

% Principal angles w.r.t R(:,1)
%R(:,13) = 