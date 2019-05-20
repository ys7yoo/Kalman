function [paramEstim LL] = em_kalman_abcd(Y, U, paramInit, max_iter, fixedParam)

if (nargin < 5)
    fixedParam = false;     % whether to fix D and R or not
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate parameter A,B,D,Q,R of the LDS:
% X(k+1) = A*X(k) + B*U(k) + v(k)
% Y(k) = C*X(k) + D*U(k) + w(k)
% Ev = 0, Evv' = Q
% Ew = 0, Eww' = R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%LL = NaN*ones(max_iter,1);
LLs = [-Inf];
converged = 0;
decreased = 0;

converged = 0;
itr = 1;
while (itr <=max_iter) && (~converged)
    %% E-step by kalman smoother
    [Xs Ps Pcs LL] = kalman_smth(Y, U, A, B, C, D, Q, R, Xo, Po);

    % check convergence 
    changeLL = (LL-LLs(end))/N;
    if ( changeLL > eps )    % check LL is increasing 
   
        %% M-step: update parameters
        % for updating A and B
        S1 = sum(Xs(:,2:end),2);
        S0 = S1+Xs(:,1)-Xs(:,end);
        SU0 = sum(U(:,1:end-1)')';

        S11 = 0; %zeros(dimX,dimX);    %Xs(:,1)*Xs(:,1)' + Ps(:,:,1);
        S10 = 0; %zeros(dimX,dimX);    %Xs(:,1)*Xo' + Pcs(:,:,1);
        S00 = 0; %zeros(dimX,dimX);    %Xo*Xo' + Po;
        SXU10 = 0;
        SXU00 = 0;
        SUU00 = 0;

        % for updating D
        SYU = Y(:,1)*U(:,1)';

        % for update R
        SYY = Y(:,1)*Y(:,1)';% + C*Ps(:,:,i)*C';
        SYX = Y(:,1)*Xs(:,1)';


        for i=2:N
            % for updating A and B
            S11 = S11 + Xs(:,i)*Xs(:,i)'+ Ps(:,:,i);
            S10 = S10 + Xs(:,i)*Xs(:,i-1)' + Pcs(:,:,i);
            S00 = S00 + Xs(:,i-1)*Xs(:,i-1)' + Ps(:,:,i-1);

            SXU10 = SXU10 + Xs(:,i)*U(:,i-1)';
            SXU00 = SXU00 + Xs(:,i-1)*U(:,i-1)';
            SUU00 = SUU00 + U(:,i-1)*U(:,i-1)';

            if (~fixedParam)
                % for updating D
                SYU = SYU + Y(:,i)*U(:,i)';


                % for update R
                SYY = SYY + Y(:,i)*Y(:,i)';% + C*Ps(:,:,i)*C';
                SYX = SYX + Y(:,i)*Xs(:,i)';
            end


    %         % altanative update rule (calculate each iteration) 
    %         % for updating Q
    %         u = Xs(:,i)-A*Xs(:,i-1)-B*D*U(:,i-1);
    %         Sx = Sx + u*u' + Ps(:,:,i) + A*Ps(:,:,i-1)*A' - Pcs(:,:,i)*A' - A*Pcs(:,:,i)';
    %         
    %         % for updating R
    %         u = Y(:,i)-C*Xs(:,i)-D*U(:,i);
    %         Sy = Sy + u*u' + C*Ps(:,:,i)*C';
        end

        % update A and B
        AB = [S10 SXU10]/[S00 SXU00; SXU00' SUU00];      %AB = S10*inv(S00);
        A = AB(:,1:dimX);
        B = AB(:,dimX+1:end);
    %     if A > 1
    %         A = .99;
    %     end
        % B = (B>1e-3).*B; % NO! B can be negative due to reversal potential

        % update D
        if (~fixedParam)
            SXU = SXU00 + Xs(:,end)*U(:,end)';
            SUU = (SUU00 + U(:,end)*U(:,end)');
            D = (SYU-C*SXU)/SUU;
            D(end) = 0;
        end


        %% update rule 1: use new AB and D for Q and R

        % update Q
        Q=(S11-AB*[S10'; SXU10'])/(N-1);   % equivalently Q=(S11-S10*inv(S00)*S10')/(N-1);
        Q=(Q+Q')/2;

        % update R
        if (~fixedParam)
            S = (S11+Xs(:,1)*Xs(:,1)'+ Ps(:,:,1));
            SXU = SXU00 + Xs(:,end)*U(:,end)';
            SUU = SUU00 + U(:,end)*U(:,end)';
            R = (SYY + C*S*C' + D*SUU*D'- D*SYU' - SYU*D' - C*SYX' -SYX*C'  ...
                 + C*SXU*D' + D*SXU'*C')/N;
            R = diag(diag(R));        % diagonalize
        end

        % update initial condition
        Xo=Xs(:,1);
        Po=Ps(:,:,1);
    
        if (R<=0)
            LL = -NaN;
            paramEstim = [];
            return;
        end
        
        
                %% check convergence here
        if ( changeLL < 1e-3) 
            disp(sprintf('Converged to %e after itr %d', LL, itr))
            converged = 1;
        end

        LLs = [LLs; LL];
        
    else % LL is decreasing ...
        disp(sprintf('[Warning] LL is decreasing from %e to %e in itr %d', LLs(end), LL, itr))
        decreased = 1;
    end 
    
    itr = itr+1;
    
end


% force non-negative!
%D = (D>0).*D;      
%D = (D>1e-3).*D;
    

paramEstim.A = A;
paramEstim.B = B;
paramEstim.C = C;
paramEstim.D = D;
paramEstim.Q = Q;
paramEstim.R = R;
paramEstim.Xo = Xo;
paramEstim.Po = Po;
paramEstim.LL = LL(end);
paramEstim.Xs = Xs;
%paramEstim.Ps = Ps;
