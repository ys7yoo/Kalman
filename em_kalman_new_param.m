
function [paramEstim paramInter] = em_kalman_new_param(Y, Iapp, ThetaInit, EM, paramTrue, saveFolderName) %#eml

PLOT_PARAMS_M_STEP = 0; % plot individual M-steps
PLOT_PARAMS_EM_ITER = 0;        % plot after all EM iterations

if (nargin < 5)
    paramTrue = [];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate parameter alpha,beta,Ke,Q,R of the LDS: with non negative
% constraint on Ke
% X(k+1) = alpha*X(k) + beta*Ke*U(k) + Vrev + v(k)
% Y(k) = C*X(k) + Ke*U(k) + w(k)
% Ev = 0, Evv' = Q
% Ew = 0, Eww' = R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = size(Y,2);
nke = EM.M; 

%UU = [U; ones(1,N)];   % already added! (May 24, 2011)
%U = UU(1:end-1,:);   % the other way. Prepare data only term
U = stackCols(Iapp,nke,0);            % July 26, 2011

%LL = NaN*ones(EM.max_iter,1);
%LLs = [-Inf];
LLs = [];

converged = 0;
decreased = 0;
reverse = 0;
itr = 1;


%% unpack initial params
alpha = ThetaInit.alpha;
beta = ThetaInit.beta;
Vrev = ThetaInit.Vrev;
Ke = ThetaInit.Ke;
Q = ThetaInit.Q;
R = ThetaInit.R;
Xo = ThetaInit.Xo;
Po = ThetaInit.Po;
C = ThetaInit.C;



%% EM iteration

if nargout>1
    % store all the intemediate params (FOR DEBUG)
    paramInter = NaN*ones(nke+5,EM.max_iter+1);
    paramInter(:,1) = [alpha; beta; Vrev; Ke'; Q; R];
    STORE_ALL_PARAM = 1;
else 
    STORE_ALL_PARAM = 0;
end

while (itr <=EM.max_iter) && ~(converged && EM.checkConverged) && ~(decreased && EM.checkDecreased) && ~reverse    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% E-step by kalman smoother
    
%     [alpha beta Vrev Ke Q R] = deal(Theta(1), Theta(2), Theta(3), Theta(4:3+nke)', Theta(4+nke), Theta(5+nke));
    
    u = Ke*U;
    %tic;
    [Xs Ps Pcs] = kalman_smth_1d(Y-Vrev, u, alpha, beta, C, 1, Q, R, Xo, Po);
    %toc; % 0.49 sec
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% M-step
    [alpha, beta, Vrev, Q, R, Ke, LL] =m_step_full_1d_new_param(Y, U, Xs, Ps, Pcs, alpha, beta, Vrev, Q, R, Ke, EM);
    %disp(sprintf('LL=%f at %dth iteration', LL, itr))
    
    if Q<0
        sprintf('[Warning] negative variance Q=%.2f',Q)
        Q = 0.1;
    end
    if R<0
        sprintf('[Warning] negative variance R=%.2f',R)
        R = 0.1;
    end
    


    if STORE_ALL_PARAM
        paramInter(:,itr+1) = [alpha; beta; Vrev; Ke'; Q; R];
    end

    
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
        fprintf('Converged to %e after itr %d: changeLL = %e\n', LL, itr, changeLL(itr))
        converged = 1;
    
    end          
   
    % increase index
    itr = itr+1;
    
end

% 
% %% plot the change of LL and params through EM iterations
% 
% if (PLOT_PARAMS_EM_ITER)
%     clf;
% 
%     % log likelihood over EM-iterations
%     subplot(331); hold on;
%     semilogy(real(LLs));
%     title (sprintf('L(\\Theta), final=%.3e',LLs(end)))
% 
% 
%     subplot(332); hold on;
%     semilogy(diff(real(LLs)));
%     if length(LLs)>1
%         title (sprintf('\\Delta L(\\Theta), final=%.3e',LLs(end)-LLs(end-1)))
%     else
%         title (sprintf('\\Delta L(\\Theta)'))
%     end
%     xlabel ('EM iteration')
% 
%     subplot(334);
%     plot(paramInter(1,1:itr), '.-');title ('\alpha');box off
%     subplot(335);
%     plot(paramInter(2,1:itr), '.-');title ('\beta');box off
%     subplot(336);
%     plot(paramInter(3,1:itr), '.-');title ('Vrev');box off
%     subplot(337);
%     plot(paramInter(4:3+nke,1:itr)', '.-');title ('Ke');box off
%     subplot(338);
%     plot(paramInter(4+nke,1:itr), '.-');title ('Q');box off
%     xlabel ('EM iteration')
%     subplot(339);
%     plot(paramInter(5+nke,1:itr), '.-');title ('R');box off
%     xlabel ('EM iteration')
% 
%     % plot true param
%     if (~isempty(paramTrue))
%         subplot(334);
%         hold on; plot(repmat(paramTrue.alpha,itr), 'r--');
%         set(gca,'ylim', [min(0.9,min(paramInter(1,1:itr))) 1])
%         subplot(335);
%         hold on; plot(repmat(paramTrue.beta,itr), 'r--');
%         set(gca,'ylim', [0 paramTrue.beta*2])
%         set(gca,'ylim', [min(paramTrue.beta*0.1,min(paramInter(2,1:itr))) max(paramTrue.beta*1.1,max(paramInter(2,1:itr)))])
%         subplot(336);
%         hold on; plot(repmat(paramTrue.Vrev,itr), 'r--');
%         %set(gca,'ylim', [-0.01 0.01])
%         set(gca,'ylim', [min(paramTrue.Vrev-0.01,min(paramInter(3,1:itr))) max(paramTrue.Vrev+0.01,max(paramInter(3,1:itr)))])
%         subplot(337);
%         hold on; plot(repmat(reshape(paramTrue.Ke,1,[]),itr,1), '--');
%         subplot(338);
%         hold on; plot(repmat(paramTrue.Q,itr), 'r--');
%         set(gca,'ylim', [0 max(paramTrue.Q*2,max(paramInter(4+nke,1:itr)))])
%         xlabel ('EM iteration')
%         subplot(339);
%         hold on; plot(repmat(paramTrue.R,itr), 'r--');
%         set(gca,'ylim', [0 max(paramTrue.R*2,max(paramInter(5+nke,1:itr)))])
%         xlabel ('EM iteration')
%     % %     subplot(259);
%     % %     hold on; plot(repmat(paramTrue.Xo,itr), 'r--');
%     % %     set(gca,'ylim', [min(paramTrue.Xo-0.1,min(paramInter(6+nke,1:itr))) max(paramTrue.Xo+0.1,max(paramInter(6+nke,1:itr)))])
%     % %     xlabel ('EM iteration')
%     % %     subplot(2,5,10);
%     % %     hold on; plot(repmat(paramTrue.Po,itr), 'r--');
%     % %     set(gca,'ylim', [0 max(paramTrue.Po*2,max(paramInter(7+nke,1:itr)))])
%     % %     xlabel ('EM iteration')    
% 
% 
%     end
% % %     % DO NOT SAVE HERE!
% % %     set(gcf, 'paperposition', [0 0 15 6]) 
% % %     set(gcf, 'papersize', [15 6]) 
% % %     saveas(1, sprintf('%s/params_EM_itr1-%d.pdf',saveFolderName,itr));
% end
% 


%% reverse parameters to previous value
if ((EM.checkDecreased && decreased) || reverse) && (length(LLs)>2)
    disp('reverse to previous params')

    LLs = LLs (1:end-1);
    LL = LLs(end);

    % reverse parameters, Xs, and LL and return
    paramEstim = paramPrev;
    paramEstim.LLs = LLs;
    paramEstim.LL = LL;

% %     % only for debugging
% %     paramEstim.Xs = XsPrev;
% %     paramEstim.Ps = PsPrev



    %keyboard
else 
    disp('Reached finial iteration')
    
    %LLs = LLs (1:end);
    
    paramEstim.alpha = alpha;
    paramEstim.beta = beta;
    paramEstim.Vrev = Vrev;
    %paramEstim.C = C;
    paramEstim.Ke = Ke;
    paramEstim.Q = Q;
    paramEstim.R = R;
    paramEstim.Xo = Xo;
    paramEstim.Po = Po;
    %paramEstim.LLes = LLes(2:end);
    paramEstim.LLs = LLs(1:end);
    paramEstim.LL = LLs(end);
% %     % only for debuggin
% %     paramEstim.Xs = Xs;
% %     paramEstim.Ps = Ps;
    %paramEstim.itr = itr-1;
    
end
 


