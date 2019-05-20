function [paramEstim] = em_kalman_abcd_new_param(Y, KeU, paramInit, EM, updateParam)

if (nargin < 4)
    % default params for EM
    EM.max_iter = 100;
    EM.checkConverged = true;
    EM.checkDecreased = true;
    EM.eps = 1e-5;
end

if (nargin < 5)    
    % determine which parameter to update
    updateParam = [1 1 0 0 1 1 1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate parameter A,B,D,Q,R of the LDS:
% X(k+1) = A*X(k) + B*KeU(k) + v(k)
% Y(k) = C*X(k) + D*KeU(k) + w(k)
% Ev = 0, Evv' = Q
% Ew = 0, Eww' = R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set initial parameters
%A = paramInit.A;
%B = paramInit.B;
C = paramInit.C;
D = paramInit.D;
Q = paramInit.Q;
R = paramInit.R;
Xo = paramInit.Xo;
Po = paramInit.Po;
Vrev = paramInit.Vrev;      % for estimating reversal potential

% parameters that really matter
aT = paramInit.A(1,:);   %only the 1st row matters 
beta = paramInit.B(1,1);
c = paramInit.C(1,1);
q = paramInit.Q(1,1);

%% check and force param conditions
%     A(2:end,1:end-1) = eye(dimX-1);
%     A(2:end,end) = 0;
%     B(2:end) = 0;
%     Q(2:end,1) = 0;
%     Q(:,2:end) = 0;
    
N = size(Y,2);
dimX = size(paramInit.A,1);

%LL = NaN*ones(max_iter,1);
LLs = -Inf;
converged = 0;
decreased = 0;

itr = 1;
%while (itr <=EM.max_iter) && (~converged)
while (itr <=EM.max_iter) && ~(converged && EM.checkConverged) && ~(decreased && EM.checkDecreased)     % for estimating reversal potential

    % build ABCD from current estimated params
    A = zeros(dimX,dimX);
    A(1,:) = aT;
    A(2:end,1:end-1) = eye(dimX-1);
    B = zeros(dimX,1);
    B(1,1) = beta;
%     %C and D is usually fixed
%     C = zeros(1,dimX);
%     C(1,1) = c;
%     
    Q = zeros(dimX,dimX);
    Q(1,1) = q;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% E-step: Kalman smoother
    [Xs Ps Pcs] = kalman_smth(Y-Vrev, KeU, A, B, C, D, Q, R, Xo, Po);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% M-step: update parameters (general version)

    % for updating A and B
    S1 = sum(Xs(:,2:end),2);
    S0 = S1+Xs(:,1)-Xs(:,end);
    SU0 = sum(KeU(:,1:end-1)')';

    P = Xs*Xs' + sum(Ps,3);
    P11 = P - (Xs(:,1)*Xs(:,1)'- Ps(:,:,1));
    PP10 = sum(Pcs(:,:,2:N),3);
    P10 = Xs(:,2:end)*Xs(:,1:end-1)' + PP10;
    PP00 = sum(Ps(:,:,1:N-1),3);  % to reuse in calculating LL
    P00 = Xs(:,1:end-1)*Xs(:,1:end-1)' + PP00;

    SXU = Xs*KeU';
    SXU00 = SXU - Xs(:,end)*KeU(:,end)';
    SXU10 = Xs(:,2:end)*KeU(:,1:end-1)';
    SUU = KeU*KeU';
    SUU00 = SUU - KeU(:,end)*KeU(:,end)';

    % for updating D and R
    SYU = (Y-Vrev)*KeU';
    SYY = (Y-Vrev)*(Y-Vrev)';
    SYX = (Y-Vrev)*Xs';

%         % update A only 
%         if (updateParam(1))
%             A = P10 /P00;      %AB = P10*inv(P00);
%         end


    % update A and B
    AB = [P10 SXU10]/[P00 SXU00; SXU00' SUU00];      %AB = P10*inv(P00);
    if (updateParam(1))
        A = AB(:,1:dimX);
    end
    if (updateParam(2))
        B = AB(:,dimX+1:end);
    end
    
    % simplified version
    ab = [P10(1,:)  SXU10(1,:)] / [P00 SXU00; SXU00' SUU00];
    aT = ab(1,1:dimX);
    beta = ab(:,dimX+1:end);
%     if A > 1
%         A = .99;
%     end
    % B = (B>1e-3).*B; % NO! B can be negative due to reversal potential

%     % update D
%     if (updateParam(4))
%         D = (SYU-C*SXU)/SUU;
%         %D(end) = 0;
%     end

    %% update rule 1: use new AB and D for Q and R

    % update Q
    % without 
    if (updateParam(5))
        Q=(P11-A*P10')/(N-1);   % equivalently Q=(P11-P10*inv(P00)*P10')/(N-1);
        % with control input 
        %Q=(P11-AB*[P10'; SXU10'])/(N-1);   % equivalently Q=(P11-P10*inv(P00)*P10')/(N-1);
        Q=(Q+Q')/2;
        
        
    end
    
    p = sum(Ps(1,1,:),3);   % sum of variances 
    p11 = p - Ps(1,1,1);    % sum of variances except the 1st term
    
    % simplified version
    x = Xs(1,:);            % xns (the 1st row vector)
    sum11 = x(1,2:end)*x(1,2:end)'+sum(Ps(1,1,2:end),3);
    sum10 = -2*aT*sum(Pcs(:,1,2:end),3);
    sum00 = aT*sum(Ps(:,:,1:end-1),3)*aT';
    q = (sum11 + sum10 + sum00  - beta^2*sum(KeU(1,2:end),2)) / (N-1);
   
    % more simplifed version using new ab 
    %q = (p11 - ab*[P10(1,:)';  SXU10(1,:)']);
    
    % update R
    if (updateParam(6))
        R = (SYY + C*P*C' + D*SUU*D'- D*SYU' - SYU*D' - C*SYX' -SYX*C'  ...
             + C*SXU*D' + D*SXU'*C')/N;
        R = diag(diag(R));        % diagonalize

%             if (R<=0)
%                 LL = -NaN;
% %                paramEstim = [];
%                 return;
%             end
    end


    %% update Vrev
    if (updateParam(7))

        Vrev = mean(Y-Xs(1,:)-KeU);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %% check convergence 
    
        %%
    % calculate LL and check convergence
% %     Theta = [alpha beta Vrev Ke Q R]';


    %% calc LL
    
    
    
    % efficiently calculate a^T Xn using filter
    axn = filter(aT(:),1,x); % test routine: clf;plot(Xs(1,:),'b');hold on; plot(axn,'g')
    temp=(x(1,2:end)-axn(1:end-1)-beta*KeU(1:end-1));
    L11 = temp*temp' + p11 -2*aT*PP10(:,1) + aT*P00*aT';
    L1 = L11 / q;
    L2 = (N-1)*log(abs(q));

    L33 = (Y-Vrev-x-KeU)*(Y-Vrev-x-KeU)' + sum(p,2);
    L3 = L33/R;
    L4 = N*log(abs(R));

% %     L5 = (Xs(1)-Xo).^2/Po;      ignore prior
% %     L6 = log(abs(Po));

    %LL = .5*(L1+L2+L3+L4+L5+L6);
    LL = -.5*(L1+L2+L3+L4);

    
    %%
    %%LL = -calcNegLLdLL_new_param (Xs, Y, U, Ps, Pcs, Theta) + .5*sum(log(squeeze(Ps)));
    disp(sprintf('LL=%f at %dth iteration', LL, itr))
    if (itr == 1)
        changeLL = Inf;
    else
        changeLL(itr) = (LL-LLs(end))/N;       
    end
    
    % store LL to LLs
    LLs = [LLs; LL];
    
    if (changeLL(itr) < -EM.eps)     % check decrease
        fprintf('[Warning] LL is decreasing from %e to %e in itr %d\n', LLs(end-1), LL, itr)
        decreased = 1;
    elseif ( changeLL(itr) < EM.eps)  % check convergence
        fprintf('Converged to %e after itr %d: changeLL = %e\n', LLs(end), itr, changeLL(itr))
        converged = 1;
    
    end 
  

    % increase index
    itr = itr+1;

end


% force non-negative!
%D = (D>0).*D;      
%D = (D>1e-3).*D;
    
% reconstruct params 
A = zeros(dimX,dimX);
A(1,:) = aT;
A(2:end,1:end-1) = eye(dimX-1);
B = zeros(dimX,1);
B(1,1) = beta;
%     %C and D is usually fixed
%     C = zeros(1,dimX);
%     C(1,1) = c;
%     
Q = zeros(dimX,dimX);
Q(1,1) = q;

paramEstim.A = A;
paramEstim.B = B;
paramEstim.C = C;
paramEstim.D = D;
paramEstim.Q = Q;
paramEstim.R = R;
paramEstim.Xo = Xo;
paramEstim.Po = Po;
paramEstim.LLs = LLs(2:end);
paramEstim.LL = LL(end);
paramEstim.Xs = Xs;
paramEstim.Ps = Ps;
paramEstim.Vrev = Vrev;





