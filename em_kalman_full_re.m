function [paramEstim LL Xs] = em_kalman_full_re(Y, U, paramInit, max_iter)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate parameter alpha,beta,Ke,Q,R of the LDS: with non negative
% constraint on Ke
% X(k+1) = alpha*X(k) + beta*Ke*U(k) + gamma + v(k)
% Y(k) = C*X(k) + Ke*U(k) + w(k)
% Ev = 0, Evv' = Q
% Ew = 0, Eww' = R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set initial parameters
alpha = paramInit.alpha;
beta = paramInit.beta;
gamma = paramInit.gamma;
C = paramInit.C;
Ke = paramInit.Ke;
Q = paramInit.Q;
R = paramInit.R;
Xo = paramInit.Xo;
Po = paramInit.Po;


N = size(Y,2);
dimX = size(alpha,1);
UU = [U; ones(1,N)];

LL = NaN*ones(max_iter,1);

converged = 0;
itr = 1;
while (itr <=max_iter) && (~converged)
    %% E-step by kalman smoother
    %[Xs Ps Pcs LL(itr)] = kalman_smth(Y, U, alpha, beta*Ke, C, Ke, Q, R, Xo, Po);
    [Xs Ps Pcs LL(itr)] = kalman_smth(Y, UU, alpha, [beta*Ke gamma], C, [Ke 0], Q, R, Xo, Po);

    
    %% M-step: update parameters
    % for updating alpha and beta
    S1 = sum(Xs(2:end));
    S0 = S1+Xs(1)-Xs(end);
    SU0 = sum(U(:,1:end-1)')';
    
    S11 = 0;    %Xs(:,1)*Xs(:,1)' + Ps(:,:,1);
    S10 = 0;    %Xs(:,1)*Xo' + Pcs(:,:,1);
    S00 = 0;    %Xo*Xo' + Po;
    SXU10 = 0;
    SXU00 = 0;
    SUU00 = 0;
   
    % for updating Ke
    SYU = Y(:,1)*U(:,1)';
    
    % for update R
    SYY = Y(:,1)*Y(:,1)';% + C*Ps(:,:,i)*C';
    SYX = Y(:,1)*Xs(:,1)';

    for i=2:N
        % for updating alpha and beta
        S11 = S11 + Xs(:,i)*Xs(:,i)'+ Ps(:,:,i);
        S10 = S10 + Xs(:,i)*Xs(:,i-1)' + Pcs(:,:,i);
        S00 = S00 + Xs(:,i-1)*Xs(:,i-1)' + Ps(:,:,i-1);
        
        SXU10 = SXU10 + Xs(:,i)*U(:,i-1)';
        SXU00 = SXU00 + Xs(:,i-1)*U(:,i-1)';
        SUU00 = SUU00 + U(:,i-1)*U(:,i-1)';

        % for updating Ke
        SYU = SYU + Y(:,i)*U(:,i)';
        
        % for update R
        SYY = SYY + Y(:,i)*Y(:,i)';% + C*Ps(:,:,i)*C';
        SYX = SYX + Y(:,i)*Xs(:,i)';

    end

    % do not update Ke and R for re-estimation
% %     %% calculate new Ke
% %     SXU = SXU00 + Xs(:,end)*U(:,end)';
% %     SUU = (SUU00 + U(:,end)*U(:,end)');
% %     %Ke = (SYU-C*SXU)/SUU;      % simple approximation
% %     Qi = inv(Q);
% %     Ri = inv(R);
% %     KeNew = (beta*Qi*(SXU10-alpha*SXU00)+Ri*(SYU-C*SXU))/(Qi*beta.^2*SUU00+Ri*SUU);
    
    
%     % update alpha and beta
%     AB = [S10 SXU10*Ke']/[S00 SXU00*Ke'; Ke*SXU00' Ke*SUU00*Ke'];      %AB = S10*inv(S00);
%     alpha = AB(:,1:dimX);
%     beta = AB(:,dimX+1:end);
    ABC = [S10 SXU10*Ke' S1]/[S00 SXU00*Ke' S0; Ke*SXU00' Ke*SUU00*Ke' Ke*SU0; S0 SU0'*Ke' N-1];
    alpha = ABC(:,1:dimX);
    beta = ABC(:,dimX+1:dimX+1);
    gamma = ABC(:,dimX+2:end);
    
    if alpha > 1
        alpha = .99;
    end

    
    %% update Q
    %Q=(S11-AB*[S10'; Ke*SXU10'])/(N-1);   % equivalently Q=(S11-S10*inv(S00)*S10')/(N-1);
     Q=(S11-ABC*[S10'; Ke*SXU10'; S1])/(N-1);   % equivalently Q=(S11-S10*inv(S00)*S10')/(N-1);
    Q=(Q+Q')/2;
    
    
% %     %% update R
% %     S = (S11+Xs(:,1)*Xs(:,1)'+ Ps(:,:,1));
% %     SXU = SXU00 + Xs(:,end)*U(:,end)';
% %     SUU = SUU00 + U(:,end)*U(:,end)';
% %     R = (SYY + C*S*C' + Ke*SUU*Ke'- Ke*SYU' - SYU*Ke' - C*SYX' -SYX*C'  ...
% %          + C*SXU*Ke' + Ke*SXU'*C')/N;
% %     R=diag(diag(R));        % diagonalize
    
    
% %     %% update Ke here
% %     Ke = (KeNew>1e-3).*KeNew;    % non-negative constraint
    
   
    %% update initial condition
    Xo=Xs(:,1);
    Po=Ps(:,:,1);
    
    if (R<=0)
        LL = -NaN;
        paramEstim = [];
        return;
    end
    itr = itr+1;
end


% force non-negative!
%Ke = (Ke>0).*Ke;      
%Ke = (Ke>1e-3).*Ke;
    

paramEstim.alpha = alpha;
paramEstim.beta = beta;
paramEstim.gamma = gamma;
paramEstim.C = C;
paramEstim.Ke = Ke;
paramEstim.Q = Q;
paramEstim.R = R;
paramEstim.Xo = Xo;
paramEstim.Po = Po;
paramEstim.LL = LL(end);
