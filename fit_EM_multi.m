function [paramEstMulti paramInit paramEstim] = fit_EM_multi(Iapp, Y, paramOneTap, FIT2)

T = length(Iapp);



clf; 
%clear paramInit paramEstMulti paramEstim LLtrace;


%% initialize params

%Need to regenrate II
M = length(paramOneTap.Ke);
II=stackRows(Iapp,M,0);

KeU=paramOneTap.Ke*II;        %% Wrong??
            
for trial = 1:FIT2.num_trial
    
    %% initialize for re-estimation with multi dim
    paramInit(trial).aT = [paramOneTap.alpha zeros(1,FIT2.dim-1)] + randn(1,FIT2.dim);
    paramInit(trial).beta = paramOneTap.beta + randn();
    paramInit(trial).c = 1;
    paramInit(trial).d = 1;
    paramInit(trial).Vrev = paramOneTap.Vrev;    
    paramInit(trial).q = 1;
    paramInit(trial).r = paramOneTap.R;
    paramInit(trial).Xo = paramOneTap.Xo;   %(paramOneTap.Xo + 2*sqrt(paramOneTap.R)*randn)*ones(FIT2.dim,1);
    paramInit(trial).Po = paramOneTap.Po*eye(FIT2.dim,FIT2.dim);   %(paramOneTap.Po + 2*sqrt(paramOneTap.R)*rand)*eye(FIT2.dim,FIT2.dim);
    

   




% %     %% Let's check smoother first
% %     A = [.99 0 0; 1 0 0; 0 1 0];
% %     B = [0.02/sum(Ke) gamma; zeros(2,2)];
% %     C = [1 0 0];
% %     D = [1 0];
% %     Q = diag ([1 0 0]);
% %     Xo = V(1)*ones(3,1);
% %     Po = 1*eye(3,3);
% %     [Xs Ps Pcs LL] = kalman_smth(Y, Uext, A, B, C, D, Q, R, Xo, Po);
% %     
% %     clf
% %     plot (V, 'r+-'); hold on
% %     plot (Xs(1,:), 'bo-')



    %%
end

%% FIT iteration
parfor trial = 1:FIT2.num_trial


    % estimate parameters
    tic
    disp(sprintf('%dth trial',trial))
    [paramEstim(trial)] = em_kalman_multi(Y(1:FIT2.N), KeU(:,1:FIT2.N), paramInit(trial), FIT2, [1 1 0 0 1 0 1]);
    disp(sprintf('  LL=%.2e (%.1fs)', paramEstim(trial).LL, toc))
% 
%     
%     subplot(121); hold on;
%     semilogy(paramEstim(trial).LLs/FIT2.N);
%     %title ('Log likelihood')
%     title ('L(\Theta)')
%     xlabel ('Iteration')
%     subplot(122); hold on;
%     plot(trial, LL(end),'x');
%     title ('Final Log likelihood')
%     xlabel ('Trial')
%     drawnow;
end
%% choose parameter set with maximum LL
[LLmax idx] = max([paramEstim(:).LL]);
paramEstMulti = paramEstim(idx);
paramEstMulti.Xs = paramEstMulti.Xs(1,:);



paramEstMulti.aT
paramEstMulti.beta
%paramEstMulti.C
%paramEstMulti.D

% 
% 
% % draw
% clf
% for trial = 1:FIT2.num_trial
%     subplot(121); hold on;
%     semilogy(real(LLtrace(trial).LL));
%     %semilogy(real([LL negLL]));
%     %title ('Log likelihood')
%     title ('L(\Theta)')
%     xlabel ('Iteration')
%     subplot(122); hold on;
%     plot(trial, LLend(trial),'x');
%     title ('Final Log likelihood')
%     xlabel ('Trial')
%     drawnow;
% end
