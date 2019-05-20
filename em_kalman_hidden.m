function [param LL] = em_kalman_hidden(Y, X, paramInit, EM, updateKinematics)

if (nargin<4)
    EM.max_iter = 10;
    EM.checkConverged = 1;
    EM.checkDecreased = 1;
end

if (nargin < 5)
    % determine which parameter to update
    updateKinematics = 1
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate parameter A,B,D,Q,R of the LDS:
% X(k+1) = A*X(k) + B*U(k) + v(k)
% Y(k) = C*X(k) + D*U(k) + w(k)
% Ev = 0, Evv' = Q
% Ew = 0, Eww' = R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set initial parameters
param.A11 = paramInit.A11;
param.A12 = paramInit.A12;
param.A21 = paramInit.A21;
param.A22 = paramInit.A22;
param.C1 = paramInit.C1;
param.C2 = paramInit.C2;
param.Q1 = paramInit.Q1;
param.Q2 = paramInit.Q2;
param.R = paramInit.R;
paramPrev = param;

% Xo = paramInit.Xo;
% Po = paramInit.Po;
param.No = paramInit.No;
param.NVo = paramInit.NVo;

N = size(Y,2)-1;
dimX = size(param.A11,1);
dimN = size(param.A21,1);

%LL = NaN*ones(max_iter,1);
LLs = [-Inf];
converged = 0;
decreased = 0;

converged = 0;
itr = 1;
%while (itr <=EM.max_iter) && (~converged)
while (itr <=EM.max_iter) && ~(converged && EM.checkConverged) && ~(decreased && EM.checkDecreased)    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% E-step: estimate hidden state only by kalman smoother
    
    YY = [X(1:dimX,2:end); Y(:,1:end-1)];      % new measurement
    UU = X(1:6,1:end-1);                    % x is now controll
    [Ns NPs NPcs LL] = kalman_smth(YY, UU, param.A22, param.A21, [param.A12; param.C2], [param.A11; param.C1], param.Q2, diag([diag(param.R); diag(param.Q1)]), param.No, param.NVo);
    param.Ns = Ns;
    param.LL = LL;
    
    
    % check convergence 
    changeLL = (LL-LLs(end))/N;   
    eps = 1e-3;
    if ( changeLL < eps && changeLL >=0) 
        disp(sprintf('Converged to %e after itr %d', LL, itr))
        converged = 1;
        LLs = [LLs; LL];
    elseif (changeLL >= eps)    % check LL is increasing
   
        %% M-step: update parameters
        % keep parameters before changing
        paramPrev = param;
        
        Xs = [X(:,1:N); Ns];

        % for updating A
        S1 = sum(Xs(2:end),2);
        S0 = S1+Xs(:,1)-Xs(:,end);

        P = Xs*Xs' + [zeros(dimX,dimX+dimN); zeros(dimN,dimX) sum(NPs,3)];
        P11 = P - Xs(:,1)*Xs(:,1)'- [zeros(dimX,dimX+dimN); zeros(dimN,dimX) NPs(:,:,1)];
        P10 = Xs(:,2:end)*Xs(:,1:end-1)' + [zeros(dimX,dimX+dimN); zeros(dimN,dimX)  sum(NPcs(:,:,2:N),3)];
        P00 = Xs(:,1:end-1)*Xs(:,1:end-1)' + [zeros(dimX,dimX+dimN); zeros(dimN,dimX) sum(NPs(:,:,1:N-1),3)];

        % for updating R
        SYY = Y(:,1:N)*Y(:,1:N)';
        SYX = Y(:,1:N)*Xs';
        
        
        % update A
        Anew = P10 /P00;      %AB = P10*inv(P00);
        if (updateKinematics)
            param.A11 = Anew(1:dimX,1:dimX);  % update kinematics or not!
        end
        param.A12 = Anew(1:dimX,dimX+1:end);
        param.A21 = Anew(dimX+1:end,1:dimX);
        param.A22 = Anew(dimX+1:end,dimX+1:end);
        
        % update C1 and C2 
        Cnew = SYX/P;
        param.C1 = Cnew(:,1:dimX);
        param.C2 = Cnew(:,dimX+1:end);
       
        % update Q
        Q=(P11-Anew*P10')/(N-1);   % equivalently Q=(P11-P10*inv(P00)*P10')/(N-1);
        %Q=(Q+Q')/2;
        Q= diag(Q);
        param.Q1 = diag(Q(1:dimX));
        param.Q2 = diag(Q(dimX+1:end));
        

        % update R
        C = [param.C1 param.C2];
        R = (SYY + C*P*C' - C*SYX' -SYX*C')/N;
        R = diag(diag(R));        % diagonalize

%             if (R<=0)
%                 LL = -NaN;
% %                paramEstim = [];
%                 return;
%             end

        % update initial condition
        No=Ns(:,1);
        
        NPo=NPs(:,:,1);
        
        LLs = [LLs; LL];
        
    else % LL is decreasing ...
        disp(sprintf('[Warning] LL is decreasing from %e to %e in itr %d', LLs(end), LL, itr))
        decreased = 1;
        
        if (EM.checkDecreased)
            % reverse parameters, Xs, and LL and return
            param = paramPrev;
            param.LLs = LLs(2:end);

            %paramEstim.Xs = kalman_smth(Y, U, paramPrev.A, paramPrev.B, paramPrev.C, paramPrev.D, paramPrev.Q, paramPrev.R, paramPrev.Xo, paramPrev.Po);
            %[Xs Ps Pcs LL] = kalman_smth(Y, U, paramPrev.A, paramPrev.B, paramPrev.C, paramPrev.D, paramPrev.Q, paramPrev.R, paramPrev.Xo, paramPrev.Po);
            %paramEstim.LL = LL; % it should be the same as LLs(end)

            LL = LLs(end);
            return
        end
        
    end 
    itr = itr+1;

end


param.LLs = LLs(2:end);