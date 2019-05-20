function [paramEstim LLs] = em_kalman(Y, Iapp, ThetaInit, EM, trial, paramTrue)

PLOT_PARAMS_M_STEP = 0; % plot individual M-steps
PLOT_PARAMS = 1;        % plot after all EM iterations

if (nargin < 5)
    trial = 0;
else
    if (PLOT_PARAMS_M_STEP)
        mkdir (sprintf('trial%d',trial));   % for DEBUG
    end
end
if (nargin < 6)
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
alpha = ThetaInit.alpha;
beta = ThetaInit.beta;
gamma = ThetaInit.gamma;
C = ThetaInit.C;
Ke = ThetaInit.Ke;
Q = ThetaInit.Q;
R = ThetaInit.R;
Xo = ThetaInit.Xo;
Po = ThetaInit.Po;


N = size(Y,2);

%UU = [U; ones(1,N)];   % already added! (May 24, 2011)
%U = UU(1:end-1,:);   % the other way. Prepare data only term
U = stackCols(Iapp,EM.M,0);            % July 26, 2011

%LL = NaN*ones(EM.max_iter,1);
%LLs = [-Inf];
LLs = [];

converged = 0;
decreased = 0;
reverse = 0;
itr = 1;

% store all the intemediate params (FOR DEBUG)
paramInter = NaN*ones(length(Ke)+5,EM.max_iter+1);
paramInter(:,1) = [alpha; beta; gamma; Ke(:); Q; R]';




while (itr <=EM.max_iter) && ~(converged && EM.checkConverged) && ~(decreased && EM.checkDecreased) && ~reverse    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% E-step by kalman smoother
    %[Xs Ps Pcs LL(itr)] = kalman_smth(Y, U, alpha, beta*Ke, C, Ke, Q, R, Xo, Po);
    
    %[Xs Ps Pcs LLe] = kalman_smth_1d(Y, UU, alpha, [beta*Ke gamma], C, [Ke 0], Q, R, Xo, Po);
    % I DON"T NEED TO PASS THIS HUGE MATRIX UU. LET"S FIX THIS
    
    % simplify
    
    %KeU = [Ke*U; ones(1,N)];
    % % KeU =  [filter(Ke,1,Iapp)'; ones(1,N)];
    %[Xs Ps Pcs LLe] = kalman_smth_1d(Y, KeU, alpha, [beta gamma], C, [1 0], Q, R, Xo, Po);
    
    % BETTER WAY TO PLUG IN GAMMA
    %KeU = [Ke*U; ones(1,N)];
    u = Ke*U;
    %tic;
    [Xs Ps Pcs LLe] = kalman_smth_1d(Y-gamma, u, alpha, beta, C, 1, Q, R, Xo, Po);
    %toc; % 0.49 sec
    %Xs = Xs + gamma;
    
    
    
% %     %% Replace with sparse implementation (NO ADVANTAGE)
% %     paramLDS.alpha = alpha;
% %     paramLDS.beta = 1/beta;
% %     paramLDS.sigx = sqrt(Q);
% %     paramLDS.sigy = sqrt(R);
% %     paramLDS.xo = 0;
% %     paramLDS.muy = gamma;
% %     
% %     tic
% %     [xhat, postVar, ~, varOff] = kalman_smth_sp(beta*u(:), Y(:), paramLDS);
% %     toc  % 0.90 sec
% %     varOff = [NaN; varOff];
% %     %
% %     clf;
% %     subplot(321);
% %     plot([Xs(:) xhat])
% %     subplot(322);
% %     plot([Xs(:)-xhat])
% %     subplot(323);
% %     plot([Ps(:) postVar])
% %     subplot(324);
% %     plot([Ps(:)-postVar])
% %     subplot(325);
% %     plot([Pcs(:) varOff])
% %     subplot(326);
% %     plot([Pcs(:)-varOff])
    
    %%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% M-step
    
        
        %keyboard
        
        %subplot(311);plot(LLs);subplot(312);plot(changeLL);subplot(313);pl
        %ot(changeLL./LLs')
    %else   %% elseif (changeLL(itr) >= eps)    % check LL is increasing
    if 1==1
        % store current params befor doing Mstep
        paramPrev.alpha = alpha;
        paramPrev.beta = beta;
        paramPrev.gamma = gamma;
        paramPrev.Ke = Ke;
        paramPrev.Q = Q;
        paramPrev.R = R;
       
% %         % only for debuggin
% %         XsPrev = Xs;
% %         PsPrev = Ps;
        

        %% initial param
        ThetaInit = [alpha; beta; gamma; Ke(:); Q; R];


        %% Update parameters according to each method 
        switch (EM.MstepMethod)
            
            case 0  % update based on partial derivatives
                
                %% changed to repeat M-step multiple times (July 28, 2011)
                ThetaNew =m_step_full_1d(Y, U, Xs, Ps, Pcs, Xo, Po, ThetaInit,  EM);
                
                % pick params at the final repetition
                ThetaNew = ThetaNew(:,end);
                
                
                


                
                %% test calling m-step multiple time (July 28, 2011)
%                 MULTIPLE_M_STEP= 5;
% % %                 params = zeros(MULTIPLE_M_STEP,22);
%                 for i=1:MULTIPLE_M_STEP
%                     tic
%                     [alpha beta gamma Ke Q R Xo Po]=m_step_full_1d(alpha, beta, gamma, Ke, Q, R, Y, U, Xo, Po, Xs, Ps, Pcs, EM.MstepConstraint);
%                     toc
% % %                     params(i,:) = [alpha beta gamma Ke Q R Xo Po]';
%                 end

% %                 box off
% %                 clf;
% %                 subplot(211);plot(params, '.--');title (sprintf('params for multiple M-step itr=%d', itr) )
% %                 subplot(212);plot(diff(params), '.--');title ('difference')
% %                 saveas(1, sprintf('multiple M-step_itr%d.pdf',itr));
                
                
                
                
            case 1  % reimplemented by Matlab opt function
                % update from opt functions without positive constraint
                Theta1 = m_step_opt(Y, U, Xs, Ps, Pcs, Xo, Po, ThetaInit, 0);

            case 2  % reimplemented by Matlab opt function with constraint mincon
                % update from opt functions without positive constraint
                Theta1 = m_step_opt(Y, U, Xs, Ps, Pcs, Xo, Po, ThetaInit, 1);
                
            case -1 % compare method (for debugging)  EM.MstepMethod = -1
                
                % update from partial derivatives 
                Theta0 =m_step_full_1d( Y, U, Xs, Ps, Pcs, Xo, Po, ThetaInit, EM);
               
                % update from opt functions without positive constraint
                ThetaInitLog = ThetaInit;
                ThetaInitLog(end-1) = -log(ThetaInitLog(end-1));
                ThetaInitLog(end) = -log(ThetaInitLog(end))
                Theta1 = m_step_opt(Y, U, Xs, Ps, Pcs, Xo, Po, ThetaInitLog, 0);
                
                % update from opt functions without positive constraint
                Theta2 = m_step_opt(Y, U, Xs, Ps, Pcs, Xo, Po, ThetaInitLog, 1);
                
                
                
                % compare LLs
                LLinit = -calcNegLLdLL (Xs, Y, U, Ps, Pcs, Xo, Po, ThetaInit);
                LL0 = -calcNegLLdLL (Xs, Y, U, Ps, Pcs, Xo, Po, Theta0(:,end));
                LL1 = -calcNegLLdLL (Xs, Y, U, Ps, Pcs, Xo, Po, Theta1);
                LL2 = -calcNegLLdLL (Xs, Y, U, Ps, Pcs, Xo, Po, Theta2);
                if(~isempty(paramTrue))
                    LL3 = -calcNegLLdLL (Xs, Y, U, Ps, Pcs, Xo, Po, [paramTrue.alpha; paramTrue.beta; paramTrue.gamma; paramTrue.Ke(:); paramTrue.Q; paramTrue.R]);
                    LLcomp = [LLinit LL0 LL1 LL2 LL3]
                else 
                    LLcomp = [LLinit LL0 LL1 LL2]
                end
                                
                
                
                % update with one of the params (Theta0)
                ThetaNew = Theta0(:,end);
                
                %ThetaNew = Theta1;
                ThetaNew = Theta2;
                ThetaNew(end-1) = exp(-ThetaNew(end-1));        % exponentiate back Q
                ThetaNew(end) = exp(-ThetaNew(end));        % exponentiate back R
%                 
                


                %% check gradient
                lossDerivFunc=@(Theta)(calcNegLLdLL (Xs, Y, U, Ps, Pcs, Xo, Po, Theta));
                lossDerivFuncLog=@(Theta)(calcNegLLparamLog (Xs, Y, U, Ps, Pcs, Xo, Po, Theta));
                
                

                % analytical gradient
                [~, dLL0] = calcNegLLdLL(Xs, Y, U, Ps, Pcs, Xo, Po, Theta0(:,end));
                [~, dLL1] = calcNegLLdLL(Xs, Y, U, Ps, Pcs, Xo, Po, Theta1(:,end));
                [~, dLL2] = calcNegLLdLL(Xs, Y, U, Ps, Pcs, Xo, Po, Theta2(:,end));


                % check gradient and Hessian 
                dx = 1e-5;
                [hess0, gr0] = compHess(lossDerivFunc, Theta0(:,end), dx);
                [hess1, gr1] = compHess(lossDerivFuncLog, Theta1, dx);
                [hess2, gr2] = compHess(lossDerivFuncLog, Theta2, dx);
                
                Cov1 = inv(hess1);
                margstd1 = sqrt(diag(Cov1));
                
                Cov2 = inv(hess2);
                margstd2 = sqrt(diag(Cov2));
                

% %                 hess2 = hess1; % hessian without B
% %                 hess2(:,4) = [];
% %                 hess2(4,:) = [];
% %                 Cov2 = inv(hess2);
% %                 margstd2 = sqrt(diag(Cov2));
% % 
% %                 %
% %                 fprintf('marginal stds:\nB = %9.2f\nke = %9.3f\nke-noB = %.3f\n', margstd1(4),mean(margstd1(6:end)),mean(margstd2(5:end)));
% % 

                
    


                                                   
                
                %%
                if (PLOT_PARAMS_M_STEP)
                    clf
                    subplot(421)
                    plot(ThetaInit(1:3), 'k'); hold on 
                    plot(Theta0(1:3,end), 'bx-'); 
                    plot(Theta1(1:3),'co-'); 
                    plot(Theta2(1:3),'ro-'); 
                    if(~isempty(paramTrue))
                        plot ([paramTrue.alpha paramTrue.beta paramTrue.gamma], 'k:')
                    end
                    title ('\alpha,  \beta, \gamma')
                    %plot(Theta0(1:3,1:end-1), ':', 'color', 0.25*[1 1 1]); 



                    subplot(422)
                    plot(ThetaInit(4:length(Ke)+3), 'k'); hold on
                    plot(Theta0(4:length(Ke)+3,end),'bx-')
                    plot(Theta1(4:length(Ke)+3), 'co-');
                    plot(Theta2(4:length(Ke)+3), 'ro-');
                    if(~isempty(paramTrue))
                        plot ([paramTrue.Ke], 'k:')
                    end
                    title ('K_e')

                    %plot(Theta0(4:length(Ke)+3,1:end-1), ':', 'color', 0.25*[1 1 1])
                    legend ('init','deriv', 'fminunc', 'fmincon', 'true')

                    subplot(423)
                    plot(ThetaInit(length(Ke)+4:length(Ke)+5),'k'); hold on
                    plot(Theta0(length(Ke)+4:length(Ke)+5,end),'bx-')
                    plot(exp(-Theta1(length(Ke)+4:length(Ke)+5)),'co-');
                    plot(exp(-Theta2(length(Ke)+4:length(Ke)+5)),'ro-');
                    if(~isempty(paramTrue))
                        plot ([paramTrue.Q paramTrue.R], 'k:')
                    end
                    title ('Q and R')


                    subplot(424)
                    plot(LLcomp, 'o-')
                    title ('LL (init, deriv, fminunc, fmincon, and true')

                    set(gcf, 'paperposition', [0 0 6 5]) 
                    set(gcf, 'papersize', [6 5]) 
                    saveas(1, sprintf('params_EM_trial%d_itr1-%d.pdf',trial,itr));


                    subplot(425); hold on
    %                 plot (dLL0(:),'b-')
    %                 plot (gr0, 'bx:'); 
                    plot (dLL1(:), 'c-'); 
                    plot (gr1, 'co:'); 
                    title ('gradient, fminunc')
                    legend ('analysis', 'numerics','location', 'NorthWest')

                    subplot(426)
                    plot (dLL2(:), 'r-');  hold on
                    plot (gr2, 'ro:'); 
                    title ('gradient, fmincon')
                    legend ('analysis', 'numerics','location', 'NorthWest')


                    % marginal stds from Hessian
                    subplot(427)
                    plot(margstd1, 'c-')

                    subplot(428)
                    plot(margstd2, 'r-')



                    set(gcf, 'paperposition', [0 0 12 15]) 
                    set(gcf, 'papersize', [12 15]) 

                    saveas(1, sprintf('trial%d/M-step_comparision_itr%d.pdf',trial,itr));    
                end

                
               
                %keyboard
                
                
        end
        
        if (Q<0 || R<0) 
             %keyboard
            sprintf('[Warning] negative variance Q=%.2f,R=%.2f',Q,R)
            reverse = 1;
        end
        
        % force non-negative!
        %Ke = (Ke>0).*Ke;      
        %Ke = (Ke>1e-3).*Ke;
        
        %% update parameter here 
        if exist('EM.fixBeta','var') & EM.fixBata
            [alpha, ~, gamma, Ke, Q, R] = deal(ThetaNew(1,end), ThetaNew(2,end), ThetaNew(3,end), ThetaNew(4:3+EM.M,end)', ...
                                             ThetaNew(4+EM.M,end),ThetaNew(5+EM.M,end));            
        else
            [alpha beta gamma Ke Q R] = deal(ThetaNew(1,end), ThetaNew(2,end), ThetaNew(3,end), ThetaNew(4:3+EM.M,end)', ...
                                             ThetaNew(4+EM.M,end),ThetaNew(5+EM.M,end));
        end
        
        % store all the intermediate params (FOR DEBUG)
        paramInter(:,itr+1) = ThetaNew;

    end  % end of M-step
    
    
    
    %%
    % calculate LL and check convergence
    Theta = [alpha beta gamma Ke Q R]';
    LL = -calcNegLLdLL (Xs, Y, U, Ps, Pcs, Xo, Po, Theta);
    %%LL = -calcNegLLdLL (Xs, Y, U, Ps, Pcs, Theta) + .5*sum(log(squeeze(Ps)));
    %disp(sprintf('LL=%f at %dth iteration', LL, itr))
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
    
    %% plot change of params during multiple M-steps
    if(PLOT_PARAMS_M_STEP && exist('params','var'))
        clf;
        %subplot(211);plot(params', '.--');title (sprintf('AllParams for multiple M-step itr=%d', itr));box off
        subplot(252);
        plot(diff(params(1,:)'), '.--');title ('\Delta \alpha');box off
        xlabel ('M-subiteration')
        subplot(253);
        plot(diff(params(2,:)'), '.--');title ('\Delta \beta');box off
        xlabel ('M-subiteration')
        subplot(254);
        plot(diff(params(3,:)'), '.--');title ('\Delta \gamma');box off
        xlabel ('M-subiteration')
        subplot(255);
        plot(diff(params(4:3+EM.M,:)'), '.--');title ('\Delta Ke');box off
        xlabel ('M-subiteration')
        subplot(257);
        plot(diff(params(4+EM.M,:)'), '.--');title ('\Delta Q');box off
        xlabel ('M-subiteration')
        subplot(258);
        plot(diff(params(5+EM.M,:)'), '.--');title ('\Delta R');box off
        xlabel ('M-subiteration')
% %         subplot(259);
% %         plot(diff(params(6+EM.M,:)'), '.--');title ('\Delta Xo');box off
% %         xlabel ('M-subiteration')
% %         subplot(2,5,10);
% %         plot(diff(params(7+EM.M,:)'), '.--');title ('\Delta Po');box off
% %         xlabel ('M-subiteration')



        % log likelihood over EM-iterations
        subplot(251); hold on;
        semilogy(real(LLs));
        title ('L(\Theta)')
        xlabel('EM iteration')

        % log likelihood during M-step
        LLduringMstep = zeros(EM.MstepRepeat+1,1);
        for i=1:EM.MstepRepeat+1
            LLduringMstep(i) = -calcNegLLdLL (Xs, Y, U, Ps, Pcs, params(:,i));
        end

        subplot(256); hold on;
        semilogy(diff(real(LLduringMstep)));
        title ('\Delta L(\Theta)')
        xlabel ('M-subiteration')


        if (mod(itr,20)==1 || itr <=20 ||converged || decreased)
            set(gcf, 'paperposition', [0 0 15 6]) 
            set(gcf, 'papersize', [15 6]) 
            saveas(1, sprintf('trial%d/multiple M-step_itr%d.pdf',trial,itr));    
        end
    end
    
    
    %% increase idx
    itr = itr+1;
    
end

% % % plot parameters during multiple M-step
% % if (converged || decreased)     % if ended by converging or decreasing LL draw again and save M-step result
% % %% plot change of params during multiple M-steps
% %     clf;
% %     %subplot(211);plot(params', '.--');title (sprintf('AllParams for multiple M-step itr=%d', itr));box off
% %     subplot(252);
% %     plot(diff(params(1,:)'), '.--');title ('\Delta \alpha');box off
% %     subplot(253);
% %     plot(diff(params(2,:)'), '.--');title ('\Delta \beta');box off
% %     subplot(254);
% %     plot(diff(params(3,:)'), '.--');title ('\Delta \gamma');box off
% %     subplot(255);
% %     plot(diff(params(4:3+EM.M,:)'), '.--');title ('\Delta Ke');box off
% %     subplot(257);
% %     plot(diff(params(4+EM.M,:)'), '.--');title ('\Delta Q');box off
% %     xlabel ('M-subiteration')
% %     subplot(258);
% %     plot(diff(params(5+EM.M,:)'), '.--');title ('\Delta R');box off
% %     xlabel ('M-subiteration')
% %     subplot(259);
% %     plot(diff(params(6+EM.M,:)'), '.--');title ('\Delta Xo');box off
% %     xlabel ('M-subiteration')
% %     subplot(2,5,10);
% %     plot(diff(params(7+EM.M,:)'), '.--');title ('\Delta Po');box off
% %     xlabel ('M-subiteration')
% % 
% % 
% %     % log likelihood over EM-iterations
% %     subplot(251); hold on;
% %     semilogy(real(LLs));
% %     title ('L(\Theta)')
% %     xlabel ('EM iteration')
% % 
% %     
% %     % log likelihood during M-step
% %     LLduringMstep = zeros(EM.MstepRepeat+1,1);
% %     for i=1:EM.MstepRepeat+1
% %         LLduringMstep(i) = -calcNegLLdLL (Xs, Y, U, Ps, Pcs, params(:,i));
% %     end
% % 
% %     subplot(256); hold on;
% %     semilogy(diff(real(LLduringMstep)));
% %     title ('\Delta L(\Theta)')
% %     xlabel ('M-subiteration')
% %                 
% %     
% %     %drawnow;
% % 
% %     set(gcf, 'paperposition', [0 0 15 6]) 
% %     set(gcf, 'papersize', [15 6]) 
% %     saveas(1, sprintf('multiple M-step_itr%d.pdf',itr));
% % end
        

%% plot param through EM iterations
if (PLOT_PARAMS)
    clf;

    % log likelihood over EM-iterations
    subplot(251); hold on;
    semilogy(real(LLs));
    title (sprintf('L(\\Theta), final=%.3e',LLs(end)))


    subplot(256); hold on;
    semilogy(diff(real(LLs)));
    if length(LLs)>1
        title (sprintf('\\Delta L(\\Theta), final=%.3e',LLs(end)-LLs(end-1)))
    else
        title (sprintf('\\Delta L(\\Theta)'))
    end
    xlabel ('EM iteration')

    subplot(252);
    plot(paramInter(1,1:itr), '.-');title ('\alpha');box off
    subplot(253);
    plot(paramInter(2,1:itr), '.-');title ('\beta');box off
    subplot(254);
    plot(paramInter(3,1:itr), '.-');title ('\gamma');box off
    subplot(255);
    plot(paramInter(4:3+EM.M,1:itr)', '.-');title ('Ke');box off
    subplot(257);
    plot(paramInter(4+EM.M,1:itr), '.-');title ('Q');box off
    xlabel ('EM iteration')
    subplot(258);
    plot(paramInter(5+EM.M,1:itr), '.-');title ('R');box off
    xlabel ('EM iteration')

    % plot true param
    if (~isempty(paramTrue))
        subplot(252);
        hold on; plot(repmat(paramTrue.alpha,itr), 'r--');
        set(gca,'ylim', [min(0.9,min(paramInter(1,1:itr))) 1])
        subplot(253);
        hold on; plot(repmat(paramTrue.beta,itr), 'r--');
        set(gca,'ylim', [0 paramTrue.beta*2])
        set(gca,'ylim', [min(paramTrue.beta*0.1,min(paramInter(2,1:itr))) max(paramTrue.beta*1.1,max(paramInter(2,1:itr)))])
        subplot(254);
        hold on; plot(repmat(paramTrue.gamma,itr), 'r--');
        %set(gca,'ylim', [-0.01 0.01])
        set(gca,'ylim', [min(paramTrue.gamma-0.01,min(paramInter(3,1:itr))) max(paramTrue.gamma+0.01,max(paramInter(3,1:itr)))])
        subplot(255);
        hold on; plot(repmat(reshape(paramTrue.Ke,1,[]),itr,1), '--');
        subplot(257);
        hold on; plot(repmat(paramTrue.Q,itr), 'r--');
        set(gca,'ylim', [0 max(paramTrue.Q*2,max(paramInter(4+EM.M,1:itr)))])
        xlabel ('EM iteration')
        subplot(258);
        hold on; plot(repmat(paramTrue.R,itr), 'r--');
        set(gca,'ylim', [0 max(paramTrue.R*2,max(paramInter(5+EM.M,1:itr)))])
        xlabel ('EM iteration')
    % %     subplot(259);
    % %     hold on; plot(repmat(paramTrue.Xo,itr), 'r--');
    % %     set(gca,'ylim', [min(paramTrue.Xo-0.1,min(paramInter(6+EM.M,1:itr))) max(paramTrue.Xo+0.1,max(paramInter(6+EM.M,1:itr)))])
    % %     xlabel ('EM iteration')
    % %     subplot(2,5,10);
    % %     hold on; plot(repmat(paramTrue.Po,itr), 'r--');
    % %     set(gca,'ylim', [0 max(paramTrue.Po*2,max(paramInter(7+EM.M,1:itr)))])
    % %     xlabel ('EM iteration')    


    end

    set(gcf, 'paperposition', [0 0 15 6]) 
    set(gcf, 'papersize', [15 6]) 
    saveas(1, sprintf('params_EM_trial%d_itr1-%d.pdf',trial,itr));
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

% %     % only for debugging
% %     paramEstim.Xs = XsPrev;
% %     paramEstim.Ps = PsPrev



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
% %     % only for debuggin
% %     paramEstim.Xs = Xs;
% %     paramEstim.Ps = Ps;
    %paramEstim.itr = itr-1;
    
end
 


