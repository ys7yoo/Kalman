
%% parameters for electrod and cell
% clear
clc

% ns = 5000;      % number of samples 
% REP = 1;


%folderName = '20120221_multi_tap_bug_fix'
%folderName = '20120308_JP_vs_YS_random_aT'
%folderName = '20120307_JP_vs_YS_same_init'
%folderName = '20120308_HH'
%folderName = '.'
folderJP = '../multiTapJP'
folderName = '.'

if ~exist(folderName,'dir')
    mkdir(folderName)
end
   


% save params

% compare with general parameterization
COMPARE_WITH_GENERAL_EM = 0;

PARAM_JP = 0;
%PARAM_JP = 1;


%% generate data using LDS
DATA_TYPE = 'GEN'
DATA_TYPE = 'HH'                %Hodgkin-Huxley neuron
%DATA_TYPE = 'LOAD'
%DATA_TYPE = 'JP'
DATA_TYPE = 'JP_DATA';          % Compare with gradient descent implemented by JP


%% Turn on/ off algorithms to run
METHOD_GRADIENT = 1;
METHOD_EM_MAT = 0;
METHOD_EM_CPP = 1;



switch DATA_TYPE
    case {'HH'}
        Ke = [0    0.4000    3.5000    1.5000    0.1500];
        
        if ~exist('ns')
            ns = input('Number of samples? ');      % number of samples 
        end
        
        I = randn(1,ns);
        
        Ke =[0 0.8 0.4 0.1 0.1 0.05];
        U = filter(Ke, 1, I);
        
        
    case {'GEN', 'LOAD'}
        
        dimX = 3;
        dimY = 1;
        dimU = 1;

        %aT = rand(1,dimX); aT = aT/sum(aT)*0.9
        rr = [.9 .7 .2]';  % specify poles (roots) of IIR filter
        avec = real(poly(rr));  % coefficients of IIR filter
        aT = -avec(2:end)
        %aT = [0.7 0.15 0.1];

        beta = 0.01;
        c = 1;
        d = 1;
        Vrev = 0;
        q = 0.01;
        r = 0.01;

        Xo = zeros(dimX,1)
        initPrior = 0.01;
        Po = initPrior*eye(dimX,dimX)
        
        if PARAM_JP
            %% use the same parameter as JP
            load (fullfile(folderJP,'multiTapJP.mat', 'avec', 'beta', 'sigx', 'sigy', 'muy', 'sigxmarg', 'dimX')) %xx yy U 

            aT = -avec(2:end)'
            dimX = length(aT);
            c=1;
            d=1;
            q = sigx^2;
            r = sigy^2;

            Xo = zeros(dimX,1);
            Po = sigxmarg^2*eye(dimX,dimX);           
            
            Vrev = muy;
        end
        
        
        

        disp ('generate data')

        


        if ~exist('ns','var')
            ns = input('number of samples? ');
        end
        if isempty(ns)
            ns = 5000;
        end

        I = randn(1,ns);
        
        Ke =[0 0.8 0.4 0.1 0.1 0.05];
        U = filter(Ke, 1, I);
        
%         U = filter([0.25 0.25 0.25 0.25], [1 0.9], I);
%         U = filter([0.25 0.25 0.25 0.25], 1, U);
%         U = filter([0.25 0.25 0.25 0.25], 1, U);
%         U = filter([0.25 0.25 0.25 0.25], 1, U);
%         U = filter([0.25 0.25 0.25 0.25], 1, U);
        %U = sin((1:ns)/25)+ 
        %U = ones(1,ns);

        A = zeros(dimX,dimX);
        A(1,:) = aT;
        A(2:end,1:end-1) = eye(dimX-1);
        B = zeros(dimX,1);
        B(1,1) = beta;
        C = [c zeros(1,dimX-1)];
        D = d;
        Q = zeros(dimX,dimX);
        Q(1,1) = q;
        R = r;


        fprintf ('generate %d samples', ns)
        [X,Y] = generate_lds(U, A, B, C, D, Q, R, Xo, Po);
        Y = Y + Vrev;
        
        
        SAVE_DATA=0
        if (SAVE_DATA)
            %% save generated data
            save('multiY.txt','Y','-ascii')
            save('multiX.txt','X','-ascii')
            save('multiU.txt','U','-ascii')
            
        end
    case 'JP_DATA'   % load date from JP's code

        if ~exist('ns')
            ns = 5000;      % number of samples 
        end
        if ~exist('REP')
            %REP = 1;
            REP = input('REP=?');
        end
        
        
        if exist('REP', 'var')
            filename = sprintf('multiTapJp_N%d_rep%d.mat',ns,REP)
        else 
            filename = sprintf('multiTapJp_N%d.mat',ns)
        end
        
        load(fullfile(folderJP,filename), 'avec', 'beta', 'sigx', 'sigy', 'xx', 'yy', 'U', 'muy', 'sigxmarg','prshat')

        % convert variables 
        aT = -avec(2:end)'
        dimX = length(aT);
        c=1;
        d=1;
        q = sigx^2
        r = sigy^2
        
        Xo = zeros(dimX,1);
        Po = sigxmarg^2*eye(dimX,dimX)
        
        Vrev = muy
        
        X = xx';   % true memtrain potential
        Y = yy';
        U = U';
        
        ns = length(Y);
        
        
        
        COMPARE_WITH_GENERAL_EM = 0;
        
        A = zeros(dimX,dimX);
        A(1,:) = aT;
        A(2:end,1:end-1) = eye(dimX-1);
        B = zeros(dimX,1);
        B(1,1) = beta;
        C = [c zeros(1,dimX-1)];
        D = d;
        Q = zeros(dimX,dimX);
        Q(1,1) = q;
        R = r;

        
end


xrange = [0 1000];
clf
if exist('I','var')
    subplot(411)
    plot(I')
    title ('input current')
    set(gca,'xlim', xrange)
end

subplot(412)
plot(U')
title ('filtered current')
set(gca,'xlim', xrange)

if exist('X','var')
    subplot(413)
    plot(X')
    title ('hidden state')
    set(gca,'xlim', xrange)
end
    %
subplot(414)
plot(Y')
title ('measurement')
set(gca,'xlim', xrange)

% save data
 



%% Kalman smoothing
% form full matrix from current params


RUN_KS = 0
if RUN_KS
    disp ('Kalman smooting')
    tic
    [Xs Ps Pcs LL] = kalman_smth(Y-Vrev, U, A, B, C, D, Q, R, Xo, Po);
    toc



    clf
    subplot(511)
    plot(U')
    title ('input')
    set(gca,'xlim', xrange)
    subplot(512)
    plot(X')
    title ('state')
    set(gca,'xlim', xrange)
    %
    subplot(513)
    plot(Y')
    title ('measurement')
    set(gca,'xlim', xrange)
    %
    subplot(514)
    plot(Xs(1,:))
    title ('estimated state')
    set(gca,'xlim', xrange)
    %
    subplot(515)
    plot(Xs(1,:)-X(1,:))
    title ('error')
    set(gca,'xlim', xrange)


    %drawnow 
    saveas(1,sprintf('KS_N%d.pdf',ns))
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% init from single tap results
% INIT_FROM_SINGLE=1
% if INIT_FROM_SINGLE
% 
% %     paramLow = [0.9*beta muy-0.1 0.9*q 0.9*r]
% %     paramHigh = [1.1*beta muy+0.1 1.1*q 1.1*r]           
%             paramLow =[0.9, 1e-6, -0.1, 0.1, 0.1];
%             paramHigh =[0.99, 0.01, 0.1, 0.5, 0.5];
% 
%     nke =6
%     [paramEst KeEst] = estimParamEM(Y, I, paramLow, paramHigh, initPrior, nke, 20, 8)
% end
% clf;
% plot ([KeEst(:) Ke(:)])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPARE ONE TRIAL FROM THE SAME INIT





COMPARE_ONE_TRIAL = 1;
% if (COMPARE_ONE_TRIAL)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's compare results of one trial for the same initialization! Feb. 16, 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate by Matlab code 
disp('estimating parameters by Matlab code')

if COMPARE_ONE_TRIAL
	EM.num_trial = 1;
else
    EM.num_trial = 32;
end

% EM.max_iter = 32;
EM.max_iter = 16;
EM.checkConverged = 1;
EM.checkDecreased = 1;
EM.dim = length(aT);
EM.eps = 1e-6;


clf; 
clear paramInit paramEst paramEstim LLtrace;

idxToEst = [1 1 0 0 1 1 1];
%idxToEst = [1 1  0 0 1 0 0];    % fix Vrev and R
idxRandomize = idxToEst
idxRandomize = [1 1  0 0 1 0 0]
idxRandomize = [1 1  1 1 1 1 0]
%idxRandomize = [0 0 0 0 0 0 0];     % start from true



%INIT_RANDOM = 0;
if ~exist('INIT_RANDOM','var')
    INIT_RANDOM = input ('Initialize randomly? (press 0 or 1)')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize 

if (idxRandomize(1))
    %paramInit.aT = aT + 0.1*randn(size(aT));
    paramInit.aT = sign(aT) + 0.1*randn(size(aT));
else
    paramInit.aT = aT;
end
if (idxRandomize(2))
    paramInit.beta = beta + 0.1*randn(size(beta)); 
else
    paramInit.beta = beta;
end
paramInit.c = c;
paramInit.d = d;
if (idxRandomize(5))
    paramInit.q = q*(1.4+0.1*randn());
else
    paramInit.q = q;
end
if (idxRandomize(6))
    paramInit.r = r*(1.4+0.1*randn());
else
    paramInit.r = r;
end

paramInit.Vrev = mean(Y); %Vrev;

paramInit.Xo = Xo;
paramInit.Po = Po; 



if METHOD_GRADIENT
    disp('Estimate parameters by gradient descent')
    
    prs0 =  [paramInit.aT(:); -log(paramInit.q); -log(paramInit.r); paramInit.Vrev; paramInit.beta];
    
    prshat = estimParamGradientDescent(Y(:), U(:), dimX, prs0)
    
end

if METHOD_EM_MAT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Now, estimate parameters
    tic    
    [paramEst] = em_kalman_multi(Y, U, paramInit, EM, idxToEst)

    % do not estimate Vrev
    %[paramEstim] = em_kalman_multi(Y, U, paramInit, EM, [1 1  0 0 1 1 0]);   

    % fix R too
    % [paramEstim] = em_kalman_multi(Y, U, paramInit, EM, [1 1  0 0 1 0 0]);   
    toc
end


if METHOD_EM_CPP

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% estimate by Cpp code 



    dim = length(aT);
    % AtLow = repmat(min(aT),1,dim)
    % AtHigh = repmat(max(aT),1,dim)
    % AtLow = -1*ones(1,dim)
    % AtHigh = 2*ones(1,dim)
    % ParamLow = [0 muy-0.1 0.01 0.01]
    % ParamHigh = [1 muy+0.1 0.1 0.1]

    % % Initialize from the  true value
    % AtLow = aT;
    % AtHigh = aT;
    % 
    % ParamLow = [beta muy q r];
    % ParamHigh = [beta muy q r];

    if (COMPARE_ONE_TRIAL)
        % start from the same init for comparison with Matlab code 
        AtLow = paramInit.aT
        AtHigh = paramInit.aT
        ParamLow = [paramInit.beta paramInit.Vrev paramInit.q paramInit.r]
        ParamHigh = [paramInit.beta paramInit.Vrev paramInit.q paramInit.r]
    else
        if INIT_RANDOM
            % random initialization
            AtLow = [0.9 -1*ones(1,dim-1)]
            AtHigh = [2*rand(1,1) 1*ones(1,dim-1)]
            ParamLow = [0 muy-0.1 0.01 0.01]
            ParamHigh = [1 muy+0.1 0.1 0.1]

        else            
            % scaled init using true
            AtLow = 0.9*aT
            AtHigh = 1.1*aT
            ParamLow = [0.9*beta muy-0.1 0.9*q 0.9*r]
            ParamHigh = [1.1*beta muy+0.1 1.1*q 1.1*r]
        end
        

        
    end

    initPrior = Po(1,1)

    tic
    [AtEstCpp, ParamEstCpp, LL] = estimParamEM_multi(Y, U, AtLow, AtHigh, ParamLow, ParamHigh, initPrior, EM.max_iter, EM.num_trial)
    toc
end

%% compare estimates 
TYPE_INIT = 'g+';
TYPE_MAT = 'bd';
TYPE_CPP='ks';
TYPE_JP='c*';

clf
row=2;
col=3;
subplot(row,col,1); hold on;
if (METHOD_EM_MAT)
    if exist('paramEstim', 'var')
        for trial = 1:length(paramEstim)
            plot(paramEstim(trial).LLs);
        end
    end
    if exist('paramEst', 'var')
        plot(paramEst.LLs,'r'); % plot the best one
    end
end

title(sprintf('LL for N=%d,itr=%d,trial=%d',ns,EM.max_iter,EM.num_trial))

subplot(row,col,2); hold on
plot (paramInit.aT(:),TYPE_INIT)
plot (aT(:), 'xr'); 
if (METHOD_EM_MAT)
    plot (paramEst.aT(:),TYPE_MAT)
end
plot (AtEstCpp,TYPE_CPP)

title ('a^T')
box off



subplot(row,col,3); hold on
plot (paramInit.beta,TYPE_INIT)
plot (beta, 'xr-'); 
if (METHOD_EM_MAT)
    plot (paramEst.beta,TYPE_MAT)
end
plot (ParamEstCpp(1),TYPE_CPP)
%set(gca, 'ylim', [0 1])
title ('\beta')
box off

legend ('init','true', 'Matlab', 'Cpp', 'Location', 'NE'); legend boxoff


subplot(row,col,4); hold on
plot (paramInit.Vrev,TYPE_INIT)
plot (Vrev, 'xr-')
if (METHOD_EM_MAT)
    plot (paramEst.Vrev,TYPE_MAT)
end
plot (ParamEstCpp(2),TYPE_CPP)
title ('V_{rev}')
%     ep = ceil(log10(abs(paramEst.Vrev)));
%set(gca, 'ylim', [-10^ep 10^ep])
% set(gca, 'ylim', [-Po(1) Po(1)])
box off

subplot(row,col,5); hold on
plot (diag(paramInit.q),TYPE_INIT)
plot (diag(q), 'xr-')
if (METHOD_EM_MAT)
    plot (diag(paramEst.q),TYPE_MAT)
end
plot (ParamEstCpp(3),TYPE_CPP)
%set(gca, 'ylim', [0 1.5*max(q,paramEst.q)])
title ('q')
box off

subplot(row,col,6); hold on
plot (diag(paramInit.r),TYPE_INIT)
plot (diag(r), 'xr-'); 
if (METHOD_EM_MAT)
    plot (diag(paramEst.r),TYPE_MAT)
end
plot (ParamEstCpp(4),TYPE_CPP)
%set(gca, 'ylim', [0 1.5*max(r,paramEst.r)])
title ('r')
box off

% if exist('REP', 'var')
%     filename = sprintf('fitting_multi_param_N%d_rep%d_T%d_cpp.pdf',ns,REP,EM.num_trial)
%     saveas(1,filename);
% else 
%     saveas(1,sprintf('fitting_multi_param_N%d_T%d_cpp.pdf',ns,EM.num_trial))
% end

% compare with JP's estimate
if (exist('prshat'))
    %% load multiTapJP.mat prshat dimX

    ahat = prshat(1:dimX);
    sigxhat = sqrt(exp(-prshat(dimX+1)));
    sigyhat = sqrt(exp(-prshat(dimX+2)));
    muyhat = prshat(dimX+3);
    betahat = prshat(dimX+4:end);


    subplot(row,col,2); hold on
    plot(ahat, TYPE_JP)
    subplot(row,col,3); hold on
    plot(betahat, TYPE_JP)
    subplot(row,col,4); hold on
    plot(muyhat, TYPE_JP)
    subplot(row,col,5); hold on
    plot(sigxhat^2, TYPE_JP)
    subplot(row,col,6); hold on
    plot(sigyhat^2, TYPE_JP)

    % replot legend
    subplot(row,col,3)
    if (METHOD_EM_MAT)
        legend ('init','true', 'Matlab', 'Cpp', 'JP', 'Location', 'NE'); legend boxoff
    else
        legend ('init','true', 'EM', 'JP', 'Location', 'NE'); legend boxoff
    end
end
subplot(row,col,3)
ylim=get(gca,'ylim');set(gca,'ylim', [0 ylim(2)]);
subplot(row,col,5)
ylim=get(gca,'ylim');set(gca,'ylim', [0 ylim(2)]);
subplot(row,col,6)
ylim=get(gca,'ylim');set(gca,'ylim', [0 ylim(2)]);




set(gcf, 'paperposition', [0 0 12 8]) 
set(gcf, 'papersize', [12 8]) 

if exist('REP', 'var')
    filename = sprintf('fitting_multi_param_N%d_rep%d_T%d_IR%d.pdf',ns,REP,EM.num_trial,INIT_RANDOM)
    saveas(1,fullfile(folderName,filename));
    filename = sprintf('fitting_multi_N%d_rep%d_T%d_IR%d.mat',ns,REP,EM.num_trial,INIT_RANDOM)
    save(fullfile(folderName,filename));
else 
    saveas(1,fullfile(folderName,sprintf('fitting_multi_param_N%d_T%d_IR%d.pdf',ns,EM.num_trial,INIT_RANDOM)))
    save(fullfile(folderName,sprintf('fitting_multi_N%d_T%d_IR%d.mat',ns,EM.num_trial,INIT_RANDOM)))
end


    

% end








%% 





% %%
% if exist('REP', 'var')
%     filename = sprintf('fitting_multi_param_N%d_rep%d_T%d.pdf',ns,REP,EM.num_trial)
%     saveas(1,filename);
% else 
%     saveas(1,sprintf('fitting_multi_param_N%d_T%d.pdf',ns,EM.num_trial))
% end


return








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate parameters by matlab code 



disp('estimating parameters...')

% EM.num_trial = 1;
EM.num_trial = 32;
EM.max_iter = 16;
EM.checkConverged = 1;
EM.checkDecreased = 1;
EM.dim = length(aT);
EM.eps = 1e-6;


clf; 
clear paramInit paramEst paramEstim LLtrace;


idxToEst = [1 1 0 0 1 1 1];
%idxToEst = [1 1  0 0 1 0 0];    % fix Vrev and R
idxRandomize = idxToEst
%idxRandomize = [1 1  0 0 1 0 0]
% initialize 
parfor trial = 1:EM.num_trial
    if (idxRandomize(1))
        paramInit(trial).aT = aT + 0.1*randn(size(aT));
    else
        paramInit(trial).aT = aT;
    end
    if (idxRandomize(2))
        paramInit(trial).beta = beta + 0.1*randn(size(beta)); 
    else
        paramInit(trial).beta = beta;
    end
    paramInit(trial).c = c;
    paramInit(trial).d = d;
    if (idxRandomize(5))
        paramInit(trial).q = q*(1+0.1*randn());
    else
        paramInit(trial).q = q;
    end
    if (idxRandomize(6))
        paramInit(trial).r = r*(1+0.1*randn());
    else
        paramInit(trial).r = r;
    end
    
    paramInit(trial).Vrev = mean(Y); %Vrev;
    
    
    paramInit(trial).Xo = Xo;
    paramInit(trial).Po = Po; 
    
    %     % start some of the value from the original value (for debug)
%     paramInit(trial).A = A;
%     paramInit(trial).B = B;
%     paramInit(trial).C = C;
%     paramInit(trial).D = D;
%     paramInit(trial).Q = Q;
%     paramInit(trial).R = R;
%     paramInit(trial).Xo = Xo;
%     paramInit(trial).Po = Po;
end
    

parfor trial = 1:EM.num_trial
    
    % trial 
    fprintf('%dth trial',trial)

       
    % estimate parameters
    tic
    
    [paramEstim(trial)] = em_kalman_multi(Y, U, paramInit(trial), EM, idxToEst);
    
    
    % do not estimate Vrev
    %[paramEstim(trial)] = em_kalman_multi(Y, U, paramInit(trial), EM, [1 1  0 0 1 1 0]);   
        
    % fix R too
    % [paramEstim(trial)] = em_kalman_multi(Y, U, paramInit(trial), EM, [1 1  0 0 1 0 0]);   
    
    toc

    
% %     %disp(sprintf('  LL=%.2e', LL(end)));
% %     subplot(121); hold on;
% %     semilogy(paramEstim(trial).LLs);
% %     %title ('Log likelihood')
% %     title ('L(\Theta)')
% %     xlabel ('Iteration')
% %     subplot(122); hold on;
% %     plot(trial, paramEstim(trial).LL,'x');
% %     title ('Final Log likelihood')
% %     xlabel ('Trial')
% %     drawnow;
end
%saveas(1, 'EM_Kalman_ABCD_trials.pdf')

% choose parameter set with maximum LL
[~, idx] = max([paramEstim(:).LL]);
paramEst = paramEstim(idx);


%% plot estimatie parameters for individual trials
numTrial = EM.num_trial;
clf
row=2;
col=3;
subplot(row,col,1); hold on;
LLinit = zeros(numTrial,1);
for trial = 1:numTrial
    LLinit(trial) = paramEstim(trial).LLs(1);
end
plot(LLinit, '+--')
plot([paramEstim(:).LL], 'o--')
plot([paramEst.LL], 'o--r')
xlabel('trial')
title ('LL')
box off

subplot(row,col,2)
%plot (aT(:), 'xr'); hold on
aTs = reshape([paramEstim(:).aT],[],EM.num_trial);
aTtrue = repmat(aT',1,EM.num_trial);
plot (aTtrue','x:'); hold on
plot (aTs','o--');
xlabel('trial')
title ('a^T')
box off


subplot(row,col,3)
plot (repmat(beta,1,numTrial), 'xr:'); hold on
plot ([paramEstim(:).beta],'o--')
set(gca, 'ylim', [0 1])
title ('\beta')
box off


subplot(row,col,4)
plot (Vrev, 'xr-'); hold on
plot (paramEst.Vrev,'o--')
title ('V_{rev}')
ep = ceil(log10(abs(paramEst.Vrev)));
%set(gca, 'ylim', [-10^ep 10^ep])
% set(gca, 'ylim', [-Po(1) Po(1)])
box off

subplot(row,col,5)
plot (repmat(q,1,numTrial), 'xr:'); hold on
plot ([paramEstim(:).q],'o--')
%set(gca, 'ylim', [0 1.5*max(q,paramEst.q)])
title ('q')
box off

subplot(row,col,6)
plot (repmat(r,1,numTrial), 'xr:'); hold on
plot ([paramEstim(:).r],'o--')
set(gca, 'ylim', [0 1.5*max(r,paramEst.r)])
title ('r')
box off


paperW = 12;
paperH = 8;
set(gcf, 'paperposition', [-0.6 0.1 paperW+1.2 paperH]) 
set(gcf, 'papersize', [paperW paperH+0.1]) 


if exist('REP', 'var')
    filename = sprintf('fitting_multi_trials_N%d_rep%d_T%d.pdf',ns,REP,EM.num_trial)
    saveas(1,filename);
else 
    saveas(1,sprintf('fitting_multi_trials_N%d_T%d.pdf',ns,EM.num_trial))
end


%% plot estimated parameters 

clf
row=2;
col=3;
subplot(row,col,1); hold on;
for trial = 1:length(paramEstim)
    plot(paramEstim(trial).LLs);
end
plot(paramEst.LLs,'r'); % plot the best one

subplot(row,col,2)
plot (aT(:), 'xr'); hold on
plot (paramEst.aT(:),'o')
% set(gca, 'ylim', [0 1])
title ('a^T')
box off


subplot(row,col,3)
plot (beta, 'xr-'); hold on
plot (paramEst.beta,'o')
%set(gca, 'ylim', [0 1])
title ('\beta')
box off


subplot(row,col,4)
plot (Vrev, 'xr-'); hold on
plot (paramEst.Vrev,'o')
title ('V_{rev}')
ep = ceil(log10(abs(paramEst.Vrev)));
%set(gca, 'ylim', [-10^ep 10^ep])
% set(gca, 'ylim', [-Po(1) Po(1)])
box off

subplot(row,col,5)
plot (diag(q), 'xr-'); hold on
plot (diag(paramEst.q),'o')
set(gca, 'ylim', [0 1.5*max(q,paramEst.q)])
title ('q')
box off

subplot(row,col,6)
plot (diag(r), 'xr-'); hold on
plot (diag(paramEst.r),'o')
set(gca, 'ylim', [0 1.5*max(r,paramEst.r)])
title ('r')
box off




% compare with JP's estimate
if (exist('prshat','var'))
    %load multiTapJP.mat prshat dimX
    
    ahat = prshat(1:dimX);
    sigxhat = sqrt(exp(-prshat(dimX+1)));
    sigyhat = sqrt(exp(-prshat(dimX+2)));
    muyhat = prshat(dimX+3);
    betahat = prshat(dimX+4:end);

    TYPE_JP='sg';
    subplot(row,col,2); hold on
    plot(ahat, TYPE_JP)
    subplot(row,col,3); hold on
    plot(betahat, TYPE_JP)
    subplot(row,col,4); hold on
    plot(muyhat, TYPE_JP)
    subplot(row,col,5); hold on
    plot(sigxhat^2, TYPE_JP)
    subplot(row,col,6); hold on
    plot(sigyhat^2, TYPE_JP)

    
end


if exist('REP', 'var')
    filename = sprintf('fitting_multi_param_N%d_rep%d_T%d.pdf',ns,REP,EM.num_trial)
    saveas(1,filename);
else 
    saveas(1,sprintf('fitting_multi_param_N%d_T%d.pdf',ns,EM.num_trial))
end




%% How about estimated voltage?
x = X(1,:);
xHat1 = paramEst.Xs(1,:);
mse1 = (x-xHat1)*(x-xHat1)' / length(x)

clf;
plot(X(1,:),'r--'); hold on
plot(paramEst.Xs(1,:),'bo-')
legend ('true', 'estim')
legend boxoff
box off
title (sprintf('mse = %.3f', mse1))
set(gca,'xlim', [0 min(length(x),200)])

paperW = 12;
paperH = 5;
set(gcf, 'paperposition', [-0.6 0.1 paperW+1.2 paperH]) 
set(gcf, 'papersize', [paperW paperH+0.1]) 




if exist('REP', 'var')
    filename = sprintf('fitting_multi_trace_N%d_rep%d_T%d.pdf',ns,REP,EM.num_trial)
    saveas(1,filename);
    filename = sprintf('fitting_multi_N%d_rep%d_T%d.mat',ns,REP,EM.num_trial)
    save(filename);
else 
    saveas(1,sprintf('fitting_multi_trace_N%d_T%d.pdf',ns,EM.num_trial))
    save(sprintf('fitting_multi_N%d_T%d.mat',ns,EM.num_trial))
end










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Re-estimate parameters with Cpp code!


% First, need to build estimParamEM_multi.mexmaci64
% See estimParamEM_multi.cpp

dim = length(aT);
% AtLow = repmat(min(aT),1,dim)
% AtHigh = repmat(max(aT),1,dim)
AtLow = -1*ones(1,dim)
AtHigh = 2*ones(1,dim)

ParamLow = [0 Vrev-0.1 0.01 0.01]
ParamHigh = [1 Vrev+0.1 0.1 0.1]


% Initialize from the  true value
AtLow = aT;
AtHigh = aT;

ParamLow = [beta Vrev q r];
ParamHigh = [beta Vrev q r];



initPrior = Po(1,1)

[AtEstCpp, ParamEstCpp] = estimParamEM_multi(Y, U, AtLow, AtHigh, ParamLow, ParamHigh, initPrior, EM.max_iter, EM.num_trial)







%% plot estimated parameters 

%load fitting_multi_N5000_T32

TYPE_CPP='s';

clf
row=2;
col=3;
subplot(row,col,1); hold on;
for trial = 1:length(paramEstim)
    plot(paramEstim(trial).LLs);
end
plot(paramEst.LLs,'r'); % plot the best one

subplot(row,col,2)
plot (aT(:), 'xr'); hold on
plot (paramEst.aT(:),'o')
plot (AtEstCpp,TYPE_CPP)

title ('a^T')
box off


subplot(row,col,3)
plot (beta, 'xr-'); hold on
plot (paramEst.beta,'o')
plot (ParamEstCpp(1),TYPE_CPP)
%set(gca, 'ylim', [0 1])
title ('\beta')
box off


subplot(row,col,4)
plot (Vrev, 'xr-'); hold on
plot (paramEst.Vrev,'o')
plot (ParamEstCpp(2),TYPE_CPP)
title ('V_{rev}')
ep = ceil(log10(abs(paramEst.Vrev)));
%set(gca, 'ylim', [-10^ep 10^ep])
% set(gca, 'ylim', [-Po(1) Po(1)])
box off

subplot(row,col,5)
plot (diag(q), 'xr-'); hold on
plot (diag(paramEst.q),'o')
plot (ParamEstCpp(3),TYPE_CPP)
%set(gca, 'ylim', [0 1.5*max(q,paramEst.q)])
title ('q')
box off

subplot(row,col,6)
plot (diag(r), 'xr-'); hold on
plot (diag(paramEst.r),'o')
plot (ParamEstCpp(4),TYPE_CPP)
%set(gca, 'ylim', [0 1.5*max(r,paramEst.r)])
title ('r')
box off

% if exist('REP', 'var')
%     filename = sprintf('fitting_multi_param_N%d_rep%d_T%d_cpp.pdf',ns,REP,EM.num_trial)
%     saveas(1,filename);
% else 
%     saveas(1,sprintf('fitting_multi_param_N%d_T%d_cpp.pdf',ns,EM.num_trial))
% end




if exist('REP', 'var')
    filename = sprintf('fitting_multi_param_N%d_rep%d_T%d_cpp.pdf',ns,REP,EM.num_trial)
    saveas(1,filename);
    filename = sprintf('fitting_multi_N%d_rep%d_T%d_cpp.mat',ns,REP,EM.num_trial)
    save(filename);
else 
    saveas(1,sprintf('fitting_multi_param_N%d_T%d_cpp.pdf',ns,EM.num_trial))
    save(sprintf('fitting_multi_N%d_T%d_cpp.mat',ns,EM.num_trial))
end
























return














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate with general method from here 

if ~COMPARE_WITH_GENERAL_EM
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate parameters with general method  and compare

disp('estimating parameters with general method')

% initialize 
for trial = 1:EM.num_trial
    
    paramInit2(trial).A = zeros(dimX,dimX);
    paramInit2(trial).A(1,:) =  paramInit(trial).aT;
    paramInit2(trial).A(2:end,1:end-1)=eye(dimX-1);
    paramInit2(trial).B = [paramInit(trial).beta; zeros(dimX-1,1)];
    paramInit2(trial).C = [paramInit(trial).c zeros(1,dimX-1)];
    paramInit2(trial).D = paramInit(trial).d;
    paramInit2(trial).Vrev = 0;
    paramInit2(trial).Q = zeros(dimX,dimX);
    paramInit2(trial).Q = paramInit(trial).q;
    paramInit2(trial).R = paramInit(trial).r;
    paramInit2(trial).Xo = paramInit(trial).Xo;
    paramInit2(trial).Po = paramInit(trial).Po;
    
    %     % start some of the value from the original value (for debug)
%     paramInit(trial).A = A;
%     paramInit(trial).B = B;
%     paramInit(trial).C = C;
%     paramInit(trial).D = D;
%     paramInit(trial).Q = Q;
%     paramInit(trial).R = R;
%     paramInit(trial).Xo = Xo;
%     paramInit(trial).Po = Po;
end
    

parfor trial = 1:EM.num_trial
    
    % trial 
    fprintf('%dth trial',trial)

       
    % estimate parameters
    tic
    
    [paramEstim2(trial)] = em_kalman_abcd(Y, U, paramInit2(trial), EM, [1 1  0 0 1 1 1]);
    
    % do not estimate Vrev
    %[paramEstim2(trial)] = em_kalman_abcd(Y, U, paramInit2(trial), EM, [1 1  0 0 1 1 0]);
    % fix R too
    % [paramEstim2(trial)] = em_kalman_abcd(Y, U, paramInit2(trial), EM, [1 1  0 0 1 0 0]);
    toc

    
% %     %disp(sprintf('  LL=%.2e', LL(end)));
% %     subplot(121); hold on;
% %     semilogy(paramEstim(trial).LLs);
% %     %title ('Log likelihood')
% %     title ('L(\Theta)')
% %     xlabel ('Iteration')
% %     subplot(122); hold on;
% %     plot(trial, paramEstim(trial).LL,'x');
% %     title ('Final Log likelihood')
% %     xlabel ('Trial')
% %     drawnow;
end
%saveas(1, 'EM_Kalman_ABCD_trials.pdf')

% choose parameter set with maximum LL
[LLmax idx] = max([paramEstim2(:).LL]);
paramEst2 = paramEstim2(idx);


%% plot estimatie parameters for individual trials
numTrial = EM.num_trial;
clf
row=2;
col=3;
subplot(row,col,1); hold on;
LLinit = zeros(numTrial,1);
for trial = 1:numTrial
    LLinit(trial) = paramEstim2(trial).LLs(1);
end
plot(LLinit, '+--')
plot([paramEstim2(:).LL], 'o--')
xlabel('trial')
title ('LL')
box off

subplot(row,col,2)
%plot (aT(:), 'xr'); hold on
As = reshape([paramEstim2(:).A],[],EM.num_trial);
aTtrue = repmat(aT',1,EM.num_trial);
plot (aTtrue','x:'); hold on
plot (As(1:3,:)','o--');
xlabel('trial')
title ('a^T')
box off


subplot(row,col,3)
plot (repmat(beta,1,numTrial), 'xr:'); hold on
Bs = [paramEstim2(:).B];
plot (Bs(1,:)','o--')
set(gca, 'ylim', [0 1])
title ('\beta')
box off


subplot(row,col,4)
plot (Vrev, 'xr-'); hold on
plot (paramEst.Vrev,'o--')
title ('V_{rev}')
ep = ceil(log10(abs(paramEst.Vrev)));
%set(gca, 'ylim', [-10^ep 10^ep])
set(gca, 'ylim', [-Po(1) Po(1)])
box off

subplot(row,col,5)
Qs = [paramEstim2(:).Q];
Qs = reshape(Qs,dimX,dimX,numTrial);
plot (repmat(q,1,numTrial), 'xr:'); hold on
plot (reshape(Qs(1,1,:),1,numTrial),'o--')
%set(gca, 'ylim', [0 1.5*max(q,paramEst.q)])
title ('q')
box off

subplot(row,col,6)
plot (repmat(r,1,numTrial), 'xr:'); hold on
plot ([paramEstim2(:).R],'o--')
set(gca, 'ylim', [0 1.5*max(r,paramEst.r)])
title ('r')
box off

saveas(1,sprintf('fitting_general_trials_N%d.pdf',ns))


%% plot estimatie two parameter sets 
clf
row=2;
col=3;

TYPE2 = 'sg'
subplot(row,col,1); hold on;
for trial = 1:length(paramEstim)
    plot(paramEstim(trial).LLs);
    plot(paramEstim2(trial).LLs,'g--');
end
plot(paramEst.LLs,'r'); % plot the best one
plot(paramEst2.LLs,'r--'); % plot the best one
box off

subplot(row,col,2)
plot (aT(:), 'xr'); hold on
plot (paramEst.aT(:),'o')
plot (paramEst2.A(1,:),TYPE2)
% set(gca, 'ylim', [-1 2])
title ('a^T')
box off

subplot(row,col,3)
plot (beta, 'xr-'); hold on
plot (paramEst.beta,'o-')
plot (paramEst2.B(:),TYPE2)
%set(gca, 'ylim', [0 1])
title ('\beta')
box off

subplot(row,col,4)
plot (Vrev, 'xr-'); hold on
plot (paramEst.Vrev,'o-')
title ('V_{rev}')
ep = ceil(log10(abs(paramEst.Vrev)));
%set(gca, 'ylim', [-10^ep 10^ep])
set(gca, 'ylim', [-Po(1) Po(1)])
box off

subplot(row,col,5)
plot (diag(q), 'xr-'); hold on
plot (diag(paramEst.q),'o-')
plot (diag(paramEst2.Q),TYPE2)
set(gca, 'ylim', [0 1.5*max(q,paramEst.q)])
title ('q')
box off

subplot(row,col,6)
plot (diag(r), 'xr-'); hold on
plot (diag(paramEst.r),'o-')
plot (diag(paramEst2.R),TYPE2)
set(gca, 'ylim', [0 1.5*max(max(r,paramEst.r),paramEst2.R)])
title ('r')
box off

saveas(1,sprintf('fitting_general_vs_multi_N%d.pdf',ns))


%% new parameterization is better! Of course. 



%% How about estimated voltage?
% x = X(1,:);
% xHat1 = paramEst.Xs(1,:);
xHat2 = paramEst2.Xs(1,:);
% mse1 = (x-xHat1)*(x-xHat1)' / length(x)
mse2 = (x-xHat2)*(x-xHat2)' / length(x)

clf;
plot(X(1,:),'r--'); hold on
plot(paramEst.Xs(1,:),'bo-')
plot(paramEst2.Xs(1,:),'gs-')
legend ('true', 'estim multi', 'estim gen')
legend boxoff
box off
title (sprintf('mse_1 = %.3f, mse_2 = %.3f', mse1, mse2))
set(gca,'xlim', [0 min(length(x),200)])


saveas(1,sprintf('fitting_general_vs_multi_trace_N%d.pdf',ns))


save(sprintf('fitting_general_vs_multi_N%d.mat',ns))


return 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [BATCH MODE] REPEAT multiple times 
%
% Copy the following command and run on melville
%
% 2012.3.1 - compare JP, EM-Matlab, EM-cpp
% 2012.1.12 - numTrial is increased to 32 and maxIter is decreased to 16
% 2012.1.5 - 

addpath Kalman


INIT_RANDOM = 0
INIT_RANDOM = 1

ns = 5000;      % number of samples 
ns = 10000;      % number of samples 
%ns = 50000;      % number of samples 

for REP = 1:50
    %%
    em_kalman_multi_test
    
    %% YS's estimate by EM
    if exist('paramEst','var')        % from Matlab
        paramEM(REP).LL = paramEst.LLs(end);
        paramEM(REP).aT = paramEst.aT;
        paramEM(REP).beta = paramEst.beta;
        paramEM(REP).Vrev = paramEst.Vrev;
        paramEM(REP).q = paramEst.q;
        paramEM(REP).r = paramEst.r;
    else                            % from Cpp
        paramEM(REP).LL = LL
        paramEM(REP).aT = AtEstCpp;
        paramEM(REP).beta = ParamEstCpp(1)
        paramEM(REP).Vrev = ParamEstCpp(2);
        paramEM(REP).q = ParamEstCpp(3);
        paramEM(REP).r = ParamEstCpp(4);
    end
    
    %% read out results from gradient descent
    if exist('prshat','var')
        %paramGD(REP).LL = LL
        paramGD(REP).aT = prshat(1:dimX);
        paramGD(REP).beta = prshat(dimX+4);
        paramGD(REP).Vrev = prshat(dimX+3);
        paramGD(REP).q = exp(-prshat(dimX+1));
        paramGD(REP).r = exp(-prshat(dimX+2));
    end
    
end

if exist('folderName','var')
    save (fullfile(folderName,sprintf('EM_vs_GD_N%d_REP%d.mat',ns, REP)), 'paramEM', 'paramGD','dimX','aT','beta','q','r','Vrev','folderName','EM','ns','INIT_RANDOM')
else
    save (sprintf('EM_vs_GD_N%d_REP%d.mat',ns, REP), 'paramEM', 'paramGD','dimX','aT','beta','q','r','Vrev','folderName','EM','ns','INIT_RANDOM')
end

% % %% [QUICK FIX]
% % for REP = 1:50  % correct parameterization
% %     paramGD(REP).q = exp(-paramGD(REP).q);
% %     paramGD(REP).r = exp(-paramGD(REP).r);
% % end
%%
    

% % %%
% % clear
% % 
% % INIT_RANDOM = 1;
% % ns = 50000;      % number of samples 
% % for REP = 1:10
% %     em_kalman_multi_test
% % end









% % %% analyis results across multiple repetitions (that are saved into
% % %% separate files)
% % 
% % clear
% % ns = 5000;
% % ns = 50000;
% % 
% % TRIAL = 32;
% % INIT_RANDOM = 0;
% % INIT_RANDOM = 1;
% % 
% % clear paramEM paramGD
% % for REP = 1:2
% %     filename = sprintf('fitting_multi_N%d_rep%d_T%d_IR%d.mat',ns,REP,TRIAL,INIT_RANDOM)
% %     %filename = 'fitting_multi_N50000_rep10_T32.mat'
% % 
% %     %load(filename,'paramEst', 'aT',' beta');
% %     load(filename);
% %     
% %     % YS's estimate by EM
% %     if exist('paramEst')        % from Matlab
% %         paramEM(REP).LL = paramEst.LLs(end);
% %         paramEM(REP).aT = paramEst.aT;
% %         paramEM(REP).beta = paramEst.beta;
% %         paramEM(REP).Vrev = paramEst.Vrev;
% %         paramEM(REP).q = paramEst.q;
% %         paramEM(REP).r = paramEst.r;
% %     else                            % from Cpp
% %         paramEM(REP).LL = LL
% %         paramEM(REP).aT = AtEstCpp;
% %         paramEM(REP).beta = ParamEstCpp(1)
% %         paramEM(REP).Vrev = ParamEstCpp(2);
% %         paramEM(REP).q = ParamEstCpp(3);
% %         paramEM(REP).r = ParamEstCpp(4);
% %         
% %     end
% %     
% %     
% %     
% %     filename = sprintf(fullfile(folderJP,'multiTapJp_N%d_rep%d.mat'),ns,REP)
% %     load(filename, 'ntaps','prshat')
% %         
% %     % JP's estimate by gradient descent
% %     paramGD(REP).ahat = prshat(1:ntaps);
% %     paramGD(REP).sigxhat =  sqrt(exp(-prshat(ntaps+1)));
% %     paramGD(REP).sigyhat = sqrt(exp(-prshat(ntaps+2)));
% %     paramGD(REP).muyhat = prshat(ntaps+3);
% %     paramGD(REP).betahat = prshat(ntaps+4:end);
% %     
% %     
% % end




% % <<<<<<< .mine  on catalix
% %     %load(filename,'paramEst', 'aT',' beta');
% %     load(filename);
% %     
% %     % YS's estimate by EM
% %     if exist('paramEst')        % from Matlab
% %         paramEM(REP).LL = paramEst.LLs(end);
% %         paramEM(REP).aT = paramEst.aT;
% %         paramEM(REP).beta = paramEst.beta;
% %         paramEM(REP).Vrev = paramEst.Vrev;
% %         paramEM(REP).q = paramEst.q;
% %         paramEM(REP).r = paramEst.r;
% %     else                            % from Cpp
% %         paramEM(REP).LL = LL
% %         paramEM(REP).aT = AtEstCpp;
% %         paramEM(REP).beta = ParamEstCpp(1)
% %         paramEM(REP).Vrev = ParamEstCpp(2);
% %         paramEM(REP).q = ParamEstCpp(3);
% %         paramEM(REP).r = ParamEstCpp(4);
% %         
% %     end
% %     
% %     
% %     
% %     filename = sprintf(fullfile(folderJP,'multiTapJp_N%d_rep%d.mat'),ns,REP)
% %     load(filename, 'dimX','prshat')
% %         
% %     % JP's estimate by gradient descent
% %     paramGD(REP).ahat = prshat(1:dimX);
% %     paramGD(REP).sigxhat =  sqrt(exp(-prshat(dimX+1)));
% %     paramGD(REP).sigyhat = sqrt(exp(-prshat(dimX+2)));
% %     paramGD(REP).muyhat = prshat(dimX+3);
% %     paramGD(REP).betahat = prshat(dimX+4:end);
% %     
% %     
% % end
% % =======
% % >>>>>>> .r1087

%% plot estimated parameters 

% combine multiple repetitions 
clear allEM allGD
allEM.aT = reshape([paramEM(:).aT],dimX,[])'
allEM.beta = [paramEM(:).beta]';
allEM.Vrev = [paramEM(:).Vrev]';
allEM.q = [paramEM(:).q]';
allEM.r = [paramEM(:).r]';

% allGD.aT = reshape([paramGD(:).ahat],dimX,[])';
% allGD.beta = [paramGD(:).betahat]';
% allGD.Vrev = [paramGD(:).muyhat]';
% allGD.q = [paramGD(:).sigxhat]';    allGD.q = allGD.q.^2;
% allGD.r = [paramGD(:).sigyhat]'; allGD.r = allGD.r.^2;

allGD.aT = reshape([paramGD(:).aT],dimX,[])';
allGD.beta = [paramGD(:).beta]';
allGD.Vrev = [paramGD(:).Vrev]';
allGD.q = [paramGD(:).q]';
allGD.r = [paramGD(:).q]';

% % <<<<<<< .mine
% % 
% % %%
% % allGD.q = exp(-allGD.q);
% % allGD.r = exp(-allGD.r);
% % 
% % %% plot estimated parameters 
% % 
% % =======
% % >>>>>>> .r1087

clf
row=2;
col=3;
cnt = 1;
% subplot(row,col,1); hold on;
% for trial = 1:length(paramEstim)
%     plot(paramEstim(trial).LLs);
% end
% plot(paramEst.LLs,'r'); % plot the best one
LINE_TYPE_TRUE = 'xr';
MARKER_SIZE_TRUE = 15;
LINE_TYPE_EM='ob';
LINE_TYPE_JP='sg';

subplot(row,col,cnt); cnt = cnt + 1;
plot (aT(:), LINE_TYPE_TRUE,'MarkerSize', MARKER_SIZE_TRUE); hold on
errorbar (mean(allEM.aT), std(allEM.aT), LINE_TYPE_EM); hold on
errorbar (mean(allGD.aT), std(allGD.aT), LINE_TYPE_JP); hold on
% set(gca, 'ylim', [0 1])
title ('a^T')
box off
legend ('true', 'EM', 'Gradient'); legend boxoff


subplot(row,col,cnt); cnt = cnt + 1;
plot (beta, LINE_TYPE_TRUE,'MarkerSize', MARKER_SIZE_TRUE); hold on
errorbar (mean(allEM.beta), std(allEM.beta), LINE_TYPE_EM); hold on
errorbar (mean(allGD.beta), std(allGD.beta), LINE_TYPE_JP); hold on
%set(gca, 'ylim', [0 1])
title ('\beta')
box off


subplot(row,col,cnt); cnt = cnt + 1;
plot (Vrev, LINE_TYPE_TRUE,'MarkerSize', MARKER_SIZE_TRUE); hold on
errorbar (mean(allEM.Vrev), std(allEM.Vrev), LINE_TYPE_EM); hold on
errorbar (mean(allGD.Vrev), std(allGD.Vrev), LINE_TYPE_JP); hold on
title ('V_{rev}')
%ep = ceil(log10(abs(paramEst.Vrev)));
%set(gca, 'ylim', [-10^ep 10^ep])
% set(gca, 'ylim', [-Po(1) Po(1)])
box off

subplot(row,col,cnt); cnt = cnt + 1;
plot (diag(q), LINE_TYPE_TRUE,'MarkerSize', MARKER_SIZE_TRUE); hold on
errorbar (mean(allEM.q), std(allEM.q), LINE_TYPE_EM); hold on
errorbar (mean(allGD.q), std(allGD.q), LINE_TYPE_JP); hold on
%set(gca, 'ylim', [0 1.5*max(q,paramEst.q)])
title ('q')
box off

subplot(row,col,cnt); cnt = cnt + 1;
plot (diag(r), LINE_TYPE_TRUE,'MarkerSize', MARKER_SIZE_TRUE); hold on
errorbar (mean(allEM.r), std(allEM.r), LINE_TYPE_EM); hold on
errorbar (mean(allGD.r), std(allGD.r), LINE_TYPE_JP); hold on
%set(gca, 'ylim', [0 1.5*max(r,paramEst.r)])
title ('r')
box off


set(gcf, 'paperposition', [0 0 12 8]) 
set(gcf, 'papersize', [12 8])
if exist ('folderName')
    saveas(1,fullfile(folderName,sprintf('EM_vs_gradient_N%d_T%d_IR%d.pdf',ns,EM.num_trial,INIT_RANDOM)))
else
    saveas(1,sprintf('EM_vs_gradient_N%d_T%d_IR%d.pdf',ns,EM.num_trial,INIT_RANDOM))
end


    
