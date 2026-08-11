function [p] = Simplex2DP(a,b,i,j)
% function [p] = Simplex2DP(a,b,i,j);
%Purpose: Evaluate 2D orthogonnal polynomial
%         on simplex at (a,b) of order (i,j).   

h1=jacobiP(i,0,0,a);
h2 = jacobiP(j,2*i+1,0,b);
p=sqrt(2.0)*h1.*h2.*(1-b).^i;
end