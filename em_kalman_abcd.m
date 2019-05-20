function [paramEstim] = em_kalman_abcd(Y, U, paramInit, EM, updateParam)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate parameter A,B,D,Q,R of the LDS:
% X(k+1) = A*X(k) + B*U(k) + v(k)
% Y(k) = C*X(k) + D*U(k) + w(k)
% Ev = 0, Evv' = Q
% Ew = 0, Eww' = R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general version (finalized on Oct. 16, 2011)

if (nargin < 4)
    % default params for EM
    EM.max_iter = 100;
    EM.checkConverged = true;
    EM.checkDecreased = true;
    EM.eps = 1e-6;
end

if (nargin < 5)    
    % determine which parameter to update
    updateParam = [1 1 0 1 1 1 1 1];
end



% set initial parameters
A = paramInit.A;
B = paramInit.B;
C = paramInit.C;
D = paramInit.D;
Q = paramInit.Q;
R = paramInit.R;
Xo = paramInit.Xo;
Po = paramInit.Po;


N = size(Y,2);
dimX = size(A,1);

%LLs = [-Inf];
LLs = [];
converged = 0;
decreased = 0;



itr = 1;
%while (itr <=EM.max_iter) && (~converged)
while (itr <=EM.max_iter) && ~(converged && EM.checkConverged) && ~(decreased && EM.checkDecreased)    

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% E-step by kalman smoother
    [Xs Ps Pcs] = kalman_smth(Y, U, A, B, C, D, Q, R, Xo, Po);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% M-step: update parameters
    % keep parameters before changing
    paramPrev.A = A;
    paramPrev.B = B;
    paramPrev.C = C;
    paramPrev.D = D;
    paramPrev.Q = Q;
    paramPrev.R = R;
    
    [A B C D Q R LL] = m_step_abcd(Y, U, Xs, Ps, Pcs, A, B, C, D, Q, R);


    
    
    
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


paramEstim.A = A;
paramEstim.B = B;
paramEstim.C = C;
paramEstim.D = D;
paramEstim.Q = Q;
paramEstim.R = R;

paramEstim.LLs = LLs;
paramEstim.LL = LLs(end);

paramEstim.Xs = Xs;
paramEstim.Ps = Ps;

% prior is just given, not to be modified
paramEstim.Xo = Xo;
paramEstim.Po = Po;

