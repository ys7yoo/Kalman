function [paramEst paramInit paramEstim Klin Vo] = fit_EM_new_param(Iapp, Y, initGuess, FIT, paramTrue, saveFolderName)
% fit EM for multiple trials
%   1. initialize params
%   2. call em_kalman_new_param.m for each trial
%
% output 
%   paramEst - best param from multiple trials
%   paramInit - param used for all the trials
%   paramEstim - param estimated for all the trials
%   Klin - estimate linear param
%   Bo - estimated DC

if nargin <5
    paramTrue = [];
end

if nargin<6
    saveFolderName = [];
else
    % make dir to save intermediate results
        
    if exist(saveFolderName,'dir')
        rmdir (saveFolderName, 's')
    end
    mkdir (saveFolderName);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate alpha,beta,Ke,Q,R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initial guess of kernel using simple linear regression
% changed to consider mean (May 24, 2011)
IIone=[stackCols(Iapp,FIT.M,0); ones(1,size(Iapp,1)) ] ;
KlinVo = IIone'\Y(:);
Klin = KlinVo(1:end-1)';
Vo = KlinVo(end);

% plot ([Klin' paramTrue.Ke'])

clear paramInit paramEstim LLtrace;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for trial = 1:FIT.num_trial
    % initialize
    if (isempty(initGuess))
        paramInit(trial).alpha = .9 + .1*rand();
        paramInit(trial).beta = .1/5.55 + 0.01*rand();
        paramInit(trial).Vrev = Vo;
        %paramInit(trial).gamma = 0.3*0.1*10.6 + 0.3*rand();
        paramInit(trial).C = 1;
        paramInit(trial).Ke = Klin + 0.1*randn(1,FIT.M);
        %paramInit(trial).Ke = [0 initGuess.sumKe/(FIT.M-2) + randn(1,FIT.M-2) 0];
    %     paramInit(trial).Q = Qinit*(.5+rand()); %min(0.25*rand+.25,Rinit);
    %     paramInit(trial).R = Rinit*(.5+rand());
        paramInit(trial).Q = initGuess.Q;
        paramInit(trial).R = initGuess.R;
        paramInit(trial).Xo = 0; 
        paramInit(trial).Po = 5;
        
        % %     % start some of the value from the original value (for debug)
% %     paramInit(trial).alpha = 0.97;
% %     paramInit(trial).beta = .1/5.55;
% %     paramInit(trial).gamma = 0.3*0.1*10.6;
% %     paramInit(trial).C = 1;    
% %     paramInit(trial).Ke =[0 0 0.4 3.5 1.5 0.15 0];
% %     paramInit(trial).Q = 1;
% %     paramInit(trial).R = 1;
% %     paramInit(trial).Xo = -60;
% %     paramInit(trial).Po = 1;
    else
        if (exist('initGuess.type'))
            switch initGuess.type
                case 'fixed'
                    paramInit(trial) = initGuess;
                    paramInit(trial).C = 1;
                %case 'center' 
            end
        else
%             %case 'range'
%             otherwise   % 
                % change initial valu with range (Mar. 13, 2011)
                if (numel(initGuess.alpha)==1)
                    paramInit(trial).alpha = min(initGuess.alpha+randn(1)*0.001,1);
                else
                    %paramInit(trial).alpha = random('unif', initGuess.alpha(1),initGuess.alpha(2));
                    paramInit(trial).alpha = initGuess.alpha(1) + rand()*(initGuess.alpha(2)-initGuess.alpha(1));
                end
                
                if FIT.fixBeta
                    paramInit(trial).beta = initGuess.beta;
                else
                    if (numel(initGuess.beta)==1)
                        paramInit(trial).beta = initGuess.beta+0.001*randn(1);
                    else
                        %paramInit(trial).beta = random('unif', initGuess.beta(1),initGuess.beta(2));
                        paramInit(trial).beta = initGuess.beta(1) + rand()*(initGuess.beta(2)-initGuess.beta(1));
                    end
                end
                
                switch numel(initGuess.Vrev)
                    case 0  % automatically set 
                        paramInit(trial).Vrev = Vo+0.001*randn(1);
                    case 1 % set fixed value
                        paramInit(trial).Vrev = initGuess.Vrev+0.001*randn(1);
                    case 2  % set a random value in the given range
                        %paramInit(trial).Vrev = random('unif', initGuess.Vrev(1),initGuess.Vrev(2));
                        paramInit(trial).Vrev = initGuess.Vrev(1) + rand()*(initGuess.Vrev(2)-initGuess.Vrev(1));
                end
                    

                paramInit(trial).C = 1;
                
                if (isempty(initGuess.Ke))
                    paramInit(trial).Ke = max(Klin + 0.1*randn(1,FIT.M),0);
                else
                    paramInit(trial).Ke = initGuess.Ke +0.01*randn(1);
                end
                if (numel(initGuess.Q)==1)      
                    paramInit(trial).Q = max(initGuess.Q + 0.01*randn(1),eps);
                else
                    %paramInit(trial).Q = random('unif', initGuess.Q(1),initGuess.Q(2));
                    paramInit(trial).Q = initGuess.Q(1) + rand()*(initGuess.Q(2)-initGuess.Q(1));
                end
                if (numel(initGuess.R)==1)      
                    paramInit(trial).R = max(initGuess.R + 0.01*randn(1),eps);
                else
                    paramInit(trial).R = initGuess.R(1) + rand()*(initGuess.R(2)-initGuess.R(1));
                end
% %                 paramInit(trial).Q = initGuess.Q;
% %                 paramInit(trial).R = initGuess.R;

                


    % %         % initGuess without noise (for debug)
    % %         paramInit(trial).alpha = initGuess.alpha;
    % %         paramInit(trial).beta = initGuess.beta;
    % %         paramInit(trial).gamma = initGuess.gamma;
    % %         paramInit(trial).Q = initGuess.Q;
    % %         paramInit(trial).R = initGuess.R;
    % %         paramInit(trial).Xo = initGuess.Xo;
    % %         paramInit(trial).Po = initGuess.Po;
    % % 
    % %         paramInit(trial).C = 1;
    % %         paramInit(trial).Ke = Klin;
        end
        
        
        % prior: There is nothing I can do for prior
        paramInit(trial).Xo = initGuess.Xo;
        paramInit(trial).Po = initGuess.Po;

           
            

    end    
end

% Now, paramInit looks like
% 
% paramInit = 
% 
% 1x8 struct array with fields:
%     alpha
%     beta
%     Vrev
%     C
%     Ke
%     Q
%     R
%     Xo
%     Po
%     



 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run FIT FIT.num_trial tiems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ticID = tic;
parfor trial = 1:FIT.num_trial
%for trial = 1:FIT.num_trial
    % estimate parameters
    
    
    tic;
    [paramEstim(trial) paramInter] = em_kalman_new_param(Y, Iapp, paramInit(trial), FIT, paramTrue, saveFolderName);        % use Iapp instead of IIone because I can easily generate Ke*Iapp with the filter command
    
    
    drawParamsItr(paramEstim(trial).LLs, paramInter);
     
    % save change of LL and params 
    set(gcf, 'paperposition', [0 0 12 10]) 
    set(gcf, 'papersize', [12 10]) 
    saveas(1, sprintf('%s/params_trial%02d.pdf',saveFolderName,trial));
   
    fprintf('%dth tiral: LL=%.2e (%d itr, %.1fs)\n', trial, paramEstim(trial).LL(end), length(paramEstim(trial)),toc)
    
    

end




%% choose parameter set with maximum LL
LLend = [paramEstim.LL];
[LLmax, idx] = max(LLend);
paramEst = paramEstim(idx);
%paramEst.Xs = XXs(idx,:);

fprintf('\nTotal %d trials: LLmax=%.2e (%.1fs)\n', FIT.num_trial, LLmax, toc(ticID))


%% now smoothe with the best param
tic;
display('filter with best param')
[paramEst.Xs paramEst.Ps] = kalman_smth_1d(Y-paramEst.Vrev, filter(paramEst.Ke,1,Iapp)', paramEst.alpha, paramEst.beta, 1, 1, paramEst.Q, paramEst.R, initGuess.Xo, initGuess.Po);
toc;


return



