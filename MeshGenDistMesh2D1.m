function [VX, VY, K, EToV,x0,y0] = MeshGenDistMesh2D1(x1,x2, y1, y2)
    % function [VX, VY, K, EToV] = MeshGenDistMesh2D()
    % Purpose : Generate 2D square mesh using DistMesh;
    % By Allan P. Engsig-Karup
    % Parameters to set/define
    % fd Distance function for mesh boundary
    % fh Weighting function for distributing elements
    % h0 Characteristic length of elements
    % Bbox Bounding box for mesh
    % param Parameters to be used in function call with DistMesh
%     fd = inline('drectangle(p,-5,5,-2,2)','p');
%     fh = @huniform;
%     h0 = 0.25;
%     Bbox = [-5 -2; 5 2];
%     PFIX = [-5 -2; 5 2; 5 -2; -5 2];
   % x1 = 0; x2 = 10;y1= 0;y2 = 10;
    fd = @(p)drectangle(p,x1,x2,y1,y2);
%     fd = inline('drectangle(p,x1,x2,y1,y2)','p');
    fh = @huniform;
    h0 = 0.25;
    Bbox = [x1 y1; x2 y2];
    PFIX = [x1 y1; x2 y2; x2 y1; x1 y2];

    x0 = (x1+x2)/2;
    y0 = (y1+y2)/2;
    param = [];
    % Call distmesh
    [Vert,EToV]=distmesh2d(fd,fh,h0,Bbox,PFIX);
    VX = Vert(:,1)'; VY = Vert(:,2)';
    Nv = length(VX); K = size(EToV,1);
    % Reorder elements to ensure counter clockwise orientation
    ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
    bx = VX(EToV(:,2)); by = VY(EToV(:,2));
    cx = VX(EToV(:,3)); cy = VY(EToV(:,3));
    D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
    i = find(D<0);
    EToV(i,:) = EToV(i,[1 3 2]);