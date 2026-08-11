function [rhsQ] = EulerRHS2D2(Q,vmapM,vmapP, vmapB, mapB,N,Dr,Ds,rx,ry,sx,sy,nx,ny, LIFT,Fscale,max,g)

    % function [rhsQ] = EulerRHS2D(Q,time, ExactSolutionBC);
    % Purpose: Evaluate RHS in 2D Euler equations, discretized on weak form
    %             with a local Lax-Friedrich flux

    %Globals2D;
    
    Nfp = N+1;
    Nfaces = 3;
    K = size(Q, 2);
    Np = ((N+1)*(N+2))/2;

    vmapM1= vmapM;
    vmapM = reshape(vmapM, Nfp*Nfaces, K); vmapP = reshape(vmapP, Nfp*Nfaces, K);

    % 1. Compute volume contributions (NOW INDEPENDENT OF SURFACE TERMS)
    [F,G,h,u,v] = EulerFluxes2D1(Q, g);

    % Compute weak derivatives
    rhsQ = zeros(size(F));
    for n=1:3
        dFdr = Dr'*F(:,:,n); dFds = Ds'*F(:,:,n);
        dGdr = Dr'*G(:,:,n); dGds = Ds'*G(:,:,n);
        rhsQ(:,:,n) = (rx.*dFdr + sx.*dFds) + (ry.*dGdr + sy.*dGds);
    end

    % 2. Compute surface contributions 
    % 2.1 evaluate '-' and '+' traces of conservative variables
    QM = zeros(Nfaces*Nfp,K,3);
    QP = zeros(Nfaces*Nfp,K,3);
    
    for n=1:3
       
        Qn = Q(:,:,n);
        QM(:,:,n) = Qn(vmapM); QP(:,:,n) = Qn(vmapP);
       
        % Set Dirichlet Boundaries conditions
       
        QM0 = reshape(QM(:,:,n), Nfaces*Nfp*K,1); % QM1 = reshape(Qn, 3*Nfp*K,1)
        QP0 = reshape(QP(:,:,n), Nfaces*Nfp*K,1);
        QP0(mapB) = - QM0(mapB);   
        QM(:,:,n) = reshape(QM0, Nfaces*Nfp,K);
        QP(:,:,n) = reshape(QP0, Nfaces*Nfp,K);
    end

    % 2.2 set boundary conditions by modifying positive traces
%     if(~isempty(ExactSolutionBC))
%       QP = feval(ExactSolutionBC, Fx, Fy, nx, ny, mapI, mapO, mapW, mapC, QP, time);
%     end

    % 2.3 evaluate primitive variables & flux functions at '-' and '+' traces
    [fM,gM,hM,uM,vM] = EulerFluxes2D1(QM, g);
    [fP,gP,hP,uP,vP] = EulerFluxes2D1(QP, g);
    
%     [fM,gM,rhoM,uM,vM,pM] = EulerFluxes2D(QM, gamma);
%     [fP,gP,rhoP,uP,vP,pP] = EulerFluxes2D(QP, gamma);

    % 2.4 Compute local Lax-Friedrichs/Rusonov numerical fluxes
%     lambda = max( sqrt(uM.^2+vM.^2) + sqrt(abs(gamma*pM./rhoM)), sqrt(uP.^2+vP.^2) + sqrt(abs(gamma*pP./rhoP)));
%     lambda = reshape(lambda, Nfp, Nfaces*K);
%     lambda = ones(Nfp, 1)*max(lambda, [], 1); 
%     lambda = reshape(lambda, Nfp*Nfaces, K);
    lambda = max*ones(Nfp*Nfaces, K);

    % 2.5 Lift fluxes
    for n=1:3
      nflux1 = nx.*(fP(:,:,n) + fM(:,:,n)) + ny.*(gP(:,:,n) + gM(:,:,n)) + lambda.*(QM(:,:,n) - QP(:,:,n));
      %nflux2 = nx.*fM(:,:,n) + ny.*gM(:,:,n);
      rhsQ(:,:,n) = rhsQ(:,:,n) + LIFT*(Fscale.*(-nflux1/2));
    end
end

