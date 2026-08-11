function [F,G,h,u,v] = EulerFluxes2D(Q, g)
     
    % function [F,G,rho,u,v,p] = EulerFluxes2D(Q, gamma)
    % Purpose: evaluate primitive variables and Euler flux functions

    % extract conserved variables
    h = Q(:,:,1);
    hu = Q(:,:,2);
    hv = Q(:,:,3);
    

    % compute primitive variables
    u = hu./h; v = hv./h;

    % compute flux functions
    F = zeros(size(Q)); 
    F(:,:,1) = hu;
    F(:,:,2) = hu.^2./h + 0.5*g*h.^2;
    F(:,:,3) = hu.*hv./h;
    
    G = zeros(size(Q));
    G(:,:,1) = hv;
    G(:,:,2) =  hu.*hv./h;
    G(:,:,3) = hv.^2./h + 0.5*g*h.^2;
    
end
