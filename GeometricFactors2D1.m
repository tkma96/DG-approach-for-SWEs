function [rx,sx,ry,sy,J, xr,yr,xs,ys] = GeometricFactors2D1(VX, VY, EToV, N)

    Np = (N+1)*(N+2)/2;

    VXT = VX(EToV);
    VYT = VY(EToV);

    xr = 0.5*(VXT(:,2) - VXT(:,1));
    yr = 0.5*(VYT(:,2) - VYT(:,1));
    xs = -0.5*(VXT(:,1) - VXT(:,3));
    ys = -0.5*(VYT(:,1) - VYT(:,3));
    J = xr.*ys - xs.*yr;
    rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;

    
    xs = repmat(xs',Np,1);
    ys = repmat(ys',Np,1);

    xr = repmat(xr',Np,1);
    yr = repmat(yr',Np,1);
    
    
    
    sx = repmat(sx',Np,1);
    sy = repmat(sy',Np,1);

    rx = repmat(rx',Np,1);
    ry = repmat(ry',Np,1);

    J = repmat(J',Np,1);

end