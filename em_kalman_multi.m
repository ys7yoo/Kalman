function [paramEstim] = em_kalman_multi(Y, KeU, paramInit, EM, updateParam)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate parameter A,B,d,Q,r of the LDS:
% X(k+1) = A*X(k) + B*KeU(k) + v(k)
% YY(k) = C*X(k) + d*KeU(k) + w(k)
% Ev = 0, Evv' = Q
% Ew = 0, Eww' = r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general version (finalized on Oct. 16, 2011)

if (nargin < 4)
    % default params for EM
    EM.max_iter = 100;
    EM.checkConverged = true;
    EM.checkDecreased = true;
    EM.eps = 1e-6;
    dimX = 3;
else
    dimX = EM.dim;
end

if (nargin < 5)    
    % determine which parameter to update
    updateParam = [1 1 0 1 1 1 1];
end



% reduced parameters
aT = paramInit.aT;
beta = paramInit.beta;
c = paramInit.c;
C = [c zeros(1,dimX-1)];
d = paramInit.d;
Vrev = paramInit.Vrev;      % for estimating reversal potential

q = paramInit.q;
r = paramInit.r;

Xo = paramInit.Xo;
Po = paramInit.Po;



N = size(Y,2);

%LLs = [-Inf];
LLs = [];
converged = 0;
decreased = 0;



itr = 1;
%while (itr <=EM.max_iter) && (~converged)
while (itr <=EM.max_iter) && ~(converged && EM.checkConverged) && ~(decreased && EM.checkDecreased)    

    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% E-step by kalman smoother
    
    % form full matrix from current params
    A = zeros(dimX,dimX);
    A(1,:) = aT;
    A(2:end,1:end-1) = eye(dimX-1);
    B = zeros(dimX,1);
    B(1,1) = beta;
    Q = zeros(dimX,dimX);
    Q(1,1) = q;
    
    YY = Y-Vrev;
    [Xs Ps Pcs] = kalman_smth(YY, KeU, A, B, C, d, Q, r, Xo, Po);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% M-step: update parameters
    % keep parameters before changing
    paramPrev.aT = aT;
    paramPrev.beta = beta;
    paramPrev.c = c;
    paramPrev.d = d;
    paramPrev.Vrev = Vrev;
    paramPrev.q = q;
    paramPrev.r = r;
      
       
    % calc square sums 
    SS00 = 0;                   
    SS10 = 0;
    for i=1:N-1
        SS00 = SS00 + Xs(:,i)*Xs(:,i)' + Ps(:,:,i); %sum(Xn*Xn') from 1 to end-1
        SS10 = SS10 + Xs(:,i+1)*Xs(:,i)' +  Pcs(:,:,i+1); %sum(X{n+1}*Xn') from 1 to end-1
    end
    %sum(Xn*Xn') from 1 to end
    SS = SS00 + (Xs(:,end)*Xs(:,end)'+ Ps(:,:,end));
    %sum(Xn*Xn') from 2 to end
    SS11 = SS - (Xs(:,1)*Xs(:,1)'+ Ps(:,:,1));   
      
    
    % calc terms that involves input
    % calc square sums 
    SXU00 = 0;                   
    SXU10 = 0;                   
    %SUU00 = 0;
    for i=1:N-1
        SXU00 = SXU00 + Xs(:,i)*KeU(:,i)'; %sum(Xn*Un') from 1 to end-1
        SXU10 = SXU10 + Xs(:,i+1)*KeU(:,i)'; %sum(X{n+1}*Un') from 1 to end-1
        %SUU00 = SUU00 + KeU(:,i)*KeU(:,i)';   %sum(Un*Un') from 1 to end-1
    end
    SXU = SXU00 + Xs(:,end)*KeU(:,end)';
    %SUU = SUU00 + KeU(:,end)*KeU(:,end)';

    % SUU** can be simplified because dim KeU = 1
    SUU00 = KeU(1,1:N-1)*KeU(1,1:N-1)';     %sum(Un*Un') from 1 to end-1
    SUU = SUU00 + KeU(1,N)*KeU(1,N)';       %sum(Un*Un') from 1 to end
    
    


    % calc terms that involves measurement (for updating d and R)
% %     SYU = 0;
% %     SYY = 0;
% %     SYX = 0;
% %     for i=1:N
% %         SYU = SYU + YY(:,i)*KeU(:,i)';
% %         SYY = SYY + YY(:,i)*YY(:,i)';
% %         SYX = SYX + YY(:,i)*Xs(:,i)';
% %     end
    
    % simplifed 
    SYU = YY(1,:)*KeU(1,:)';
    SYY = YY(1,:)*YY(1,:)';
    SYX = YY(1,:)*Xs(1,:)';
    
    
    
%     % calc terms that involves measurement (for updating d and r)
%     syu = YY(1,:)*KeU(1,:)';
%     syy = YY(1,:)*YY(1,:)';
%     syx = YY(1,:)*Xs(1,:)';
    

    % calc quadratic terms
    quad1 = SS11(1,1)+aT*SS00*aT'+beta*SUU00*beta'  - 2*aT*SS10(1,:)' + 2*aT*SXU00*beta' -2*SXU10(1,:)*beta' ;
    %quad2 = SYY + C*SS*C' + d*SUU*d'- d*SYU' - SYU*d' - 2*c*SYX' + 2*C*SXU*d';
    % simplify for C=[1 0 0 ...] d=1
    quad2 = SYY + SS(1,1) + SUU - 2*SYU - 2*SYX' + 2*SXU(1,1);
         
         

    % update A and B
    AB = [SS10(1,:) SXU10(1,:)]/[SS00 SXU00; SXU00' SUU00];
    if ~sum(isnan(AB(:)))
        if (updateParam(1))
            aT = AB(:,1:dimX);
        end
        if (updateParam(2))
            beta = AB(:,dimX+1:end);
        end
    end    

%     % update d
%     if (updateParam(4))
%         d = (syu-C*sxu)/suu;
%     end

    % update Q    
    if (updateParam(5))
        q = quad1 / (N-1);
    end

    % update r
    if (updateParam(6))
        r = (quad2)/N;
    end

    %% update Vrev
    if (updateParam(7))

        Vrev = mean(Y-c*Xs(1,:)-d*KeU);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% calc LL & check convergence 
    L1 =  trace(quad1)/q + (N-1)*log(abs(q));
    L2 =  trace(quad2)/r + N*log(abs(det(r)));
    
    LL = -0.5*(L1 + L2  );
% %     %disp(sprintf('LL=%f at %dth iteration', LL, itr))



    
    
     

    % check convergence 
    if (itr == 1)
        changeLL = Inf;
    else
        changeLL = (LL-LLs(itr-1))/N;       
    end
    
    if (changeLL < -EM.eps)     % check decrease
        fprintf('[Warning] LL is decreasing from %e to %e at itr %d\n', LLs(end), LL, itr)
        decreased = 1;
        
        if (EM.checkDecreased)
            disp('reversing params')
            % reverse parameters, Xs, and LL and return
            paramEstim = paramPrev;
            
            paramEstim.LLs = LLs;
            paramEstim.LL = LLs(end);
            
            paramEstim.Xs = Xs;
            paramEstim.Ps = Ps;
            
            % prior is just given, not to be modified
            paramEstim.Xo = Xo;
            paramEstim.Po = Po;

            return
        end
        
    elseif ( changeLL < EM.eps)  % check convergence
        fprintf('Converged to %e at itr %d: changeLL = %e\n', LL, itr, changeLL)
        converged = 1;
    end
    
    % store LL to LLs
    LLs = [LLs; LL];

    
    itr = itr+1;

end


paramEstim.aT = aT;
paramEstim.beta = beta;
paramEstim.c = c;
paramEstim.d = d;
paramEstim.Vrev = Vrev;
paramEstim.q = q;
paramEstim.r = r;

paramEstim.LLs = LLs;
paramEstim.LL = LLs(end);

paramEstim.Xs = Xs;
paramEstim.Ps = Ps;

% prior is just given, not to be modified
paramEstim.Xo = Xo;
paramEstim.Po = Po;
