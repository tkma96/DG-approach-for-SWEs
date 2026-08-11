% Purpose : Setup script, building operators, grid, metric, and connectivity tables.
% Definition of constants
N = 3;
Nfp = N+1; Np = (N+1)*(N+2)/2; Nfaces=3; NODETOL = 1e-12;

% Compute nodal set
[r,s] = Nodes2D(N);% [r,s] = xytors(x,y);

% Build reference element matrices
V = Vandermonde2D(N,r,s);
invV = inv(V);
MassMatrix = invV'*invV;
[Dr,Ds, Vr, Vs] = Dmatrices2D(N,r,s,V);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('MS_MatricesPrecomputation.mat')
V = MS_MatricesPrecomputation(N).V2D;
Dr= MS_MatricesPrecomputation(N).Da;
Ds = MS_MatricesPrecomputation(N).Db;
invM = V*V';
Drw = invM*Dr'/invM;
Dsw = invM*Ds'/invM;

%%%%%%%%%%%%%%%MeshgenDishMesh2D%%%%%%%%%%%%
% x1 = -5; x2 = 5; y1 = -5; y2 = 5;
% [VX, VY, K, EToV,x0,y0] = MeshGenDistMesh2D1(x1,x2, y1, y2);
[VX, VY, K, EToV] = MeshGenDistMesh2D();

% build coordinates of all the nodes
va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
x = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));
y = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));

% Set up the initial condition
x0 = 0; y0 = 0;
[Q0] = intialsolution(x0,y0 ,x,y);

%plot the initial h

for i=1:size(x,2)
    figure(1)
    trig = delaunay(x(:,i), y(:,i));
    trisurf(trig, x(:,i), y(:,i), Q0(:,i,1))
    hold on
end
shading interp


% find all the nodes that lie on each edge
fmask1   = find( abs(s+1) < NODETOL)'; 
fmask2   = find( abs(r+s) < NODETOL)';
fmask3   = find( abs(r+1) < NODETOL)';
Fmask  = [fmask1;fmask2;fmask3]';
Fx = x(Fmask(:), :); Fy = y(Fmask(:), :);
% 
% % Create surface integral terms
 LIFT = Lift2D(Fmask,s,r,N,V);

% calculate geometric factors
[rx,sx,ry,sy,J, xr,
    
yr,xs,ys] = GeometricFactors2D1(VX, VY, EToV, N);

% calculate geometric factors
[nx, ny, sJ] = Normals2D1(xr,yr,xs,ys,Fmask,N);
Fscale = sJ./(J(Fmask,:));44

% Build connectivity matrix
[EToE, EToF,BCType] = tiConnect2D(EToV);

% Build connectivity maps
[mapM, mapP, vmapM, vmapP, vmapB, mapB] = BuildMaps2D(EToV, EToE, EToF, VX, VY, x, y, Fmask, N, NODETOL);
% % 
% Build boundaries conditons %%

In = 1; Out = 2; Wall = 3; Far = 4; Cyl = 5; Dirichlet = 6;
Neuman = 7; Slip = 8;
BCType = Dirichlet*BCType;
[mapI,mapO,mapW,mapF,mapC,mapD,mapN,mapS,vmapI,vmapO,vmapW,vmapF,vmapC,vmapD,vmapN,vmapS] = BuildBCMaps2D(vmapM,BCType, N,In, Out, Wall, Far, Cyl,Dirichlet, Neuman,Slip);


max = 0;
g=9.8;
% [rhsQ] = EulerRHS2D1(Q0,vmapM,vmapP, vmapB, mapB,N,Dr,Ds,rx,ry,sx,sy,nx,ny, LIFT,Fscale,max,g);
% % % Compute weak operators (could be done in preprocessing to save time)
%Drw = (V*Vr')/(V*V'); Dsw = (V*Vs')/(V*V');

ht = 0.01;
Qn = Q0;
for i=1:2
    K1 = EulerRHS2D1(Qn,vmapM,vmapP, vmapB, mapB,N,Dr,Ds,rx,ry,sx,sy,nx,ny, LIFT,Fscale,max,g);
    K2 = EulerRHS2D1(Qn+K1*ht/2,vmapM,vmapP, vmapB, mapB,N,Dr,Ds,rx,ry,sx,sy,nx,ny, LIFT,Fscale,max,g);
    K3 = EulerRHS2D1(Qn+K2*ht/2,vmapM,vmapP, vmapB, mapB,N,Dr,Ds,rx,ry,sx,sy,nx,ny, LIFT,Fscale,max,g);
    K4 = EulerRHS2D1(Qn+K3*ht,vmapM,vmapP, vmapB, mapB,N,Dr,Ds,rx,ry,sx,sy,nx,ny, LIFT,Fscale,max,g);
%     [rhsQ] = EulerRHS2D2(Q,vmapM,vmapP, vmapB, mapB,N,Dr,Ds,rx,ry,sx,sy,nx,ny, LIFT,Fscale,max,g);
    Qn = Qn + ht*(K1+2*K2+2*K3+K4)/6;
end
%plot final solution
for i=1:size(x,2)
    figure(1)
    trig = delaunay(x(:,i), y(:,i));
    trisurf(trig, x(:,i), y(:,i), Q(:,i,1))
    hold on
end
shading interp
