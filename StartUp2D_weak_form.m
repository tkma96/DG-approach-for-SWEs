N = 4;
%[r,s] = Nodes2D(N);% [r,s] = xytors(x,y);
zz = ((N+1)*(N+2))/2;
% Build reference element matrices
load('MS_MatricesPrecomputation.mat')
V = MS_MatricesPrecomputation(N).V2D;
Dr= MS_MatricesPrecomputation(N).Da;
Ds = MS_MatricesPrecomputation(N).Db;
V1 = MS_MatricesPrecomputation(N).V1D;
invM = V*V';
Drw = invM*Dr'/invM;
Dsw = invM*Ds'/invM;
r = MS_MatricesPrecomputation(N).a;
s = MS_MatricesPrecomputation(N).b;
for i=1:zz
   
trisurf(r(:,i), s(:,i), V(:,1))
    
end 

% line([-1,1],[-1,-1],'color', 'black','LineWidth',3);
% hold on
% line([-1,-1], [-1, 1],'color', 'black','LineWidth',3);
% line([1,-1],[-1,1],'color', 'black','LineWidth',3);
% plot(r,s, '.','MarkerFaceColor','red','MarkerEdgeColor','red', 'Markersize',90);

