function [Q] = intialsolution(x0,y0 ,x,y)
A = 1;
sig_X = 0.25;
sig_Y = 0.2;
Np = size(x,1);
K = size(x,2);
Q = zeros(Np, K, 3);
f = @(X,Y) A*exp(-(X-x0).^2/(2*sig_X^2) - ((Y-y0).^2/(2*sig_Y^2)));

Q(:,:,1) = f(x,y);
    
end
    
    
    
    
    