function [paramEstim LLs] = em_kalman_1d(Y, U, paramInit, EM, paramTrue)

if (nargin < 5)
    paramTrue = [];
end

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

UU = [U; ones(1,N)];

%LL = NaN*ones(EM.max_iter,1);
%LLs = [-Inf];
LLs = [];

converged = 0;
decreased = 0;
reverse = 0;
itr = 1;
while (itr <=EM.max_iter) && ~(converged && EM.checkConverged) && ~(decreased && EM.checkDecreased) && ~reverse
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% E-step by kalman smoother
    %[Xs Ps Pcs LL(itr)] = kalman_smth(Y, U, alpha, beta*Ke, C, Ke, Q, R, Xo, Po);
    [Xs Ps Pcs LLe] = kalman_smth_1d(Y, UU, alpha, [beta*Ke gamma], C, [Ke 0], Q, R, Xo, Po);
    %disp(sprintf('LLe=%f',LLe))
    
    
    % calculate LL => moved to after M-step (Jul. 26, 2010)
    Theta = [alpha beta gamma Ke Q R Xo Po]';
    LL = -calcNegLLdLL (Xs, Y, U, Ps, Pcs, Theta) + .5*sum(log(squeeze(Ps)));              
    if (itr == 1)
        changeLL = Inf;
    else
        changeLL(itr) = (LL-LLs(end))/N;       
    end
    % store LL to LLs
    LLs = [LLs; LL];
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% M-step
    eps = 1e-6;
    if ( changeLL(itr) < eps && changeLL(itr) >= 0)  % check convergence
        disp(sprintf('Converged to %e after itr %d', LLs(end), itr))
        LL
        changeLL
        changeLL./LLs'
        
        converged = 1;
    elseif (changeLL(itr)<0)     % check decrease
        disp(sprintf('[Warning] LL is decreasing from %e to %e in itr %d', LLs(end-1), LL, itr))
        LL
        changeLL
        changeLL./LLs'
        
        decreased = 1;
        
        %keyboard
        
        %subplot(311);plot(LLs);subplot(312);plot(changeLL);subplot(313);pl
        %ot(changeLL./LLs')
    else   %% elseif (changeLL(itr) >= eps)    % check LL is increasing
        % store current params befor doing Mstep
        paramPrev.alpha = alpha;
        paramPrev.beta = beta;
        paramPrev.gamma = gamma;
        paramPrev.Ke = Ke;
        paramPrev.Q = Q;
        paramPrev.R = R;
        paramPrev.Xo = Xo;
        paramPrev.Po = Po;
        XsPrev = Xs;
        PsPrev = Ps;
        
        switch (EM.MstepMethod)
            case 0  % update based on partial derivatives
                [alpha beta gamma Ke Q R Xo Po]=m_step_full_1d(alpha, beta, gamma, Ke, Q, R, Y, U, Xo, Po, Xs, Ps, Pcs, EM.MstepConstraint);
                
            case 1  % reimplemented by Matlab opt function
                %[alpha beta gamma Ke Q R Xo Po]=m_step_opt(alpha, beta, gamma, Ke, Q, R, Y, U, Xo, Po, Xs, Ps, Pcs, EM.MstepConstraint);
                [alpha beta gamma Ke Q R Xo Po]=m_step_opt(alpha, beta, gamma, Ke, Q, R, Y, U, Xo, Po, Xs, Ps, Pcs, 0);

            case 2  % reimplemented by Matlab opt function with constraint mincon
                %[alpha beta gamma Ke Q R Xo Po]=m_step_opt(alpha, beta, gamma, Ke, Q, R, Y, U, Xo, Po, Xs, Ps, Pcs, EM.MstepConstraint);
                [alpha beta gamma Ke Q R Xo Po]=m_step_opt(alpha, beta, gamma, Ke, Q, R, Y, U, Xo, Po, Xs, Ps, Pcs, 1);
%                 Theta = [alpha beta gamma Ke Q R Xo Po]';
%                 LL = -negLL;
            case -1 % compare method (for debugging)  EM.MstepMethod = -1
                
                % update from partial derivatives 
                %[alpha beta gamma Ke Q R Xo Po]=m_step_full(alpha, beta, gamma, Ke, Q, R, Y, U, Xo, Po, Xs, Ps, Pcs, EM.MstepConstraint);
                [alpha beta gamma Ke Q R Xo Po]=m_step_full(alpha, beta, gamma, Ke, Q, R, Y, U, Xo, Po, Xs, Ps, Pcs, 0);
                Theta0 = [alpha beta gamma Ke Q R Xo Po]';
                LL0 = -calcNegLLdLL (Xs, Y, U, Ps, Pcs, Theta0);
                
                % update from opt functions
                %[alpha beta gamma Ke Q R Xo Po negLL]=m_step_opt(alpha, beta, gamma, Ke, Q, R, Y, U, Xo, Po, Xs, Ps, Pcs, EM.MstepConstraint);
                [alpha beta gamma Ke Q R Xo Po]=m_step_opt(alpha, beta, gamma, Ke, Q, R, Y, U, Xo, Po, Xs, Ps, Pcs, 0);
                Theta1 = [alpha beta gamma Ke Q R Xo Po]'
                LL1 = -calcNegLLdLL (Xs, Y, U, Ps, Pcs, Theta1);
                                
                [LLe LL LL0 LL1]
                
                %%
                clf
                subplot(221)
                plot([Theta0(1:3) Theta1(1:3)],'x-')
                subplot(222)
                plot([Theta0(4:length(Ke)+3) Theta1(4:length(Ke)+3)],'x-')
                legend ('deriv', 'opt')
                subplot(223)
                plot([Theta0(length(Ke)+4:length(Ke)+5) Theta1(length(Ke)+4:length(Ke)+5)],'x-')
                subplot(224)
                %plot([Theta0(length(Ke)+6) Theta1(length(Ke)+6)],'x-')
                plot([Theta0(length(Ke)+6:end) Theta1(length(Ke)+6:end)],'x-')
                
                keyboard
                
                
        end
        
        if (Q<0 || R<0) 
             %keyboard
            sprintf('[Warning] negative variance Q=%.2f,R=%.2f',Q,R)
            reverse = 1;
        end
        
        % force non-negative!
        %Ke = (Ke>0).*Ke;      
        %Ke = (Ke>1e-3).*Ke;
        

        
     end
    
    itr = itr+1;
end

%% reverse parameters to previous value
if ((EM.checkDecreased && decreased) || reverse) && (length(LLs)>2)
    disp('reverse to previous params')

    LLs = LLs (1:end-1);
    LL = LLs(end);

    % reverse parameters, Xs, and LL and return
    paramEstim = paramPrev;
    paramEstim.LLs = LLs;
    paramEstim.LL = LL;

    paramEstim.Xs = XsPrev;
    paramEstim.Ps = PsPrev



    %keyboard
else 
    disp('Reached finial iteration')
    
    %LLs = LLs (1:end);
    
    paramEstim.alpha = alpha;
    paramEstim.beta = beta;
    paramEstim.gamma = gamma;
    %paramEstim.C = C;
    paramEstim.Ke = Ke;
    paramEstim.Q = Q;
    paramEstim.R = R;
    paramEstim.Xo = Xo;
    paramEstim.Po = Po;
    %paramEstim.LLes = LLes(2:end);
    paramEstim.LLs = LLs(1:end);
    paramEstim.LL = LLs(end);
    paramEstim.Xs = Xs;
    paramEstim.Ps = Ps;
    %paramEstim.itr = itr-1;
    
end
 


