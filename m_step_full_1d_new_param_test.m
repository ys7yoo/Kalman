%% typical parameters 
% from 
% Dataset:Cell_04_15_2010_BD_n0.1
% Frame = 1 - 8

clear
%alpha0 = 0.998; 
alpha0 = 0.99; 
%beta0  = 0.017; 
beta0  = 0.01; 
%Vrev0 = -67.593;
Vrev0 = 0;
% Q0     = 1e-04;
% R0     = 5e-02;
Q0     = 1e-01;
R0     = 1e-01;
Ke0 = [ -0.0029    3.0602    3.3801    2.2578    1.9500    1.7128    1.4858    1.3227    1.1470    1.0245    0.9104    0.8011    0.7300    0.6347 ...
    0.5797    0.5126    0.4712    0.4176    0.3597    0.3365    0.3158    0.2751    0.2438    0.2208    0.2189    0.1891    0.1547    0.1549 ...
    0.1462    0.1309    0.1169    0.0968    0.1095    0.0861    0.0841    0.0701    0.0617    0.0732    0.0585    0.0550    0.0443    0.0396 ...
    0.0419    0.0348    0.0326    0.0178    0.0147    0.0428    0.0159    0.0096]
    

%Ke0 = Ke0(1:5:end); % for shorter
Ke0 = [0 0.8 0.4 0.2 0.1 0.05];


Xo = 0;
Po = 0.01;

%% generate state and measurement
addpath ../YsTools
ns = 10000;
Iapp = 10*randn(1,ns);
% Iapp = filter([0.25 0.25 0.25 0.25], [1 0.9], Iapp);
% Iapp = filter([0.25 0.25 0.25 0.25], 1, Iapp);
% Iapp = filter([0.25 0.25 0.25 0.25], 1, Iapp);
% Iapp = filter([0.25 0.25 0.25 0.25], 1, Iapp);
% Iapp = filter([0.25 0.25 0.25 0.25], 1, Iapp);

nke = length(Ke0);
II = stackCols(Iapp,nke,0);            % July 26, 2011
U = Ke0*II;

%%

[X,Y] = generate_lds(U, alpha0, beta0, 1, 1, Q0, R0, Xo, Po);
Y = Y+Vrev0;
clf
subplot(311)
plot(X)
subplot(312)
plot(Y)


%% KS with true param
[Xs Ps Pcs] = kalman_smth_1d(Y-Vrev0, Ke0*II, alpha0, beta0, 1, 1, Q0, R0, Xo, Po);

mse = (X-Xs)*(X-Xs)'/ns  

clf;
subplot(211)
plot(Y); title ('measurement')
subplot(212)
plot([X' Xs'])
title (sprintf('V_m (mse=%.3f)',mse))

%% Or, read from previous data 


% load X.txt
% load Y.txt
% 
% %load II.txt
% load Iapp.txt
% nke = 6
% II = stackCols(Iapp,nke,0);            % July 26, 2011

% to save 
Xt = X';
Yt = Y';
It = Iapp';
save ('X.txt', '-ascii','Xt')
save ('Y.txt', '-ascii', 'Yt')
save ('I.txt', '-ascii', 'It')


    
%% 0. initialize 
clc
EM.num_trial = 1;
EM.M = nke;
% EM.num_trial = 8;
EM.max_iter = 50;
% EM.max_iter = 1;
EM.checkConverged = 1;
EM.checkDecreased = 1;
EM.eps = 1e-6;
%
EM.fixBeta = 1;
EM.MstepRepeat = 1;
EM.MstepConstraint = 0;


% % % true param
% % paramTrue = [alpha0; beta0; Vrev0; Q0; R0];
% % 
% % 
% % % initialize randomly
% % alpha = 0.9+ 0.1*rand();
% % beta = beta0*2*rand();
% % % beta = 0.1*rand();
% % 
% % % alpha = alpha0; 
% % % beta = beta0; 
% % 
% % 
% % % Q = Q0; 
% % Q = rand();
% % % R = R0;
% % R = rand();
 
 
% initialize from one param
alpha = 0.8
beta = 0.01
Vrev = 0
Q = 0.5
R = 0.5





% Ke = Ke0;
% Vrev = Vrev0;  
% initialize Ke and Vrev from linear 
IIone=[II; ones(1,size(II,2))];
KlinVo = IIone'\Y(:);
Ke = KlinVo(1:end-1)'
Vrev = KlinVo(end);

clf
plot([Ke0' Ke']);   % compare with true Ke
legend('true', 'linear')

% store initial param
params = [alpha; beta; Vrev; Q; R];
Kes = Ke(:);
LLs = NaN;
mses=NaN;

%
% EM.max_iter = 1;
for itr = 1:EM.max_iter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% E-step 
    % Call E-step from either Matlab version or Cpp version
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% 1. Kalman smoothe with wrong param (params)
    %[alpha beta Vrev Q R] = deal(params(1,end), params(2,end), params(3,end), params(4,end), params(5,end));
    [Xs Ps Pcs] = kalman_smth_1d(Y-Vrev, Ke*II, alpha, beta, 1, 1, Q, R, Xo, Po);

    if exist('X','var')
        mses(itr+1) = (X-Xs)*(X-Xs)'/ns;
    end

%     clf;
%     plot([X' Xs'])

    %% 2. M-step 

%     % save to compare with cpp implementation 
%     save('m_step.mat', 'Y', 'II', 'Xs', 'Ps', 'Pcs', 'alpha', 'beta', 'Vrev', 'Q', 'R', 'Ke');
%     save -ascii Y.txt Y
%     save -ascii II.txt II
%     save -ascii Xs.txt Xs
%     save -ascii Ps.txt Ps
%     save -ascii Pcs.txt Pcs
%     save -ascii alpha.txt alpha
%     save -ascii beta.txt beta
%     save -ascii Vrev.txt Vrev
%     save -ascii Q.txt Q
%     save -ascii R.txt R
%     save -ascii Ke.txt Ke
%     params = [alpha, beta, Vrev, Q, R,];
%     save -ascii params.txt params
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% M-step 
    % Call M-step from either Matlab version or Cpp version
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matlab version
%     [alpha, beta, Vrev, Q, R, Ke, LL] = m_step_full_1d_new_param(Y, II, Xs, Ps, Pcs, alpha, beta, Vrev, Q, R, Ke, EM);
    
    % call cpp version
    [alpha, beta, Vrev, Q, R, Ke, LL] = m_step_1d(Y, Iapp, Xs, Ps, Pcs, alpha, beta, Vrev, Q, R, Ke, EM);
    
    
    
    % store params
    params(1,itr+1) = alpha;
    params(2,itr+1) = beta;
    params(3,itr+1) = Vrev;
    params(4,itr+1) = Q;
    params(5,itr+1) = R;
    Kes(:,itr+1) = Ke(:);
    LLs(itr+1) = LL;
    
    
    clf;
    subplot(511)
    plot(Y); title ('measurement')
    subplot(512)
   if exist('X','var')
        plot([X' Xs'])
        title (sprintf('V_m (mse=%.3f)',mses(itr+1)))
   else
       plot([Xs'])
       title (sprintf('V_m (mse=%.3f)'))
   end
    subplot(525)
    if exist('paramTrue','var')
        plot([paramTrue(1:5,end) params(1:5,end)], '.-')
    else
        plot(params(1:5,end), '.-')
    end
    subplot(526) 
    plot(params')% change of params over iteration
    subplot(527)
    plot([Ke0(:) Kes(:,end)],'.-')
    subplot(528)
    plot(Kes')
    subplot(5,2,9)
    plot(LLs)
    subplot(5,2,10)
    plot(mses)
    


    drawnow
end

% plot change of params
alpha
beta
Vrev
Q
R
Ke


%% KS with estimated param
% [Xs Ps Pcs] = kalman_smth_1d(Y-Vrev, Ke*II, alpha, beta, 1, 1, Q, R, Xo, Po);
% mse = (X-Xs)*(X-Xs)'/ns  

[Xp Pp Xf Pf] = kalman_filt_1d(Y-Vrev, Ke*II, alpha, beta, 1, 1, Q, R, Xo, Po);
mse = (X-Xf)*(X-Xf)'/ns  



clf;
subplot(211)
plot(Y); title ('measurement')
subplot(212)
plot([X' Xs'])
title (sprintf('V_m (mse=%.3f)',mse))


return


%% compare above code with cpp implementation
% Input: Y, Iapp, initial params

load Xs.txt
load Ps.txt
load Pcs.txt

clf
subplot(311)
plot([X(:) Xs(:)])
subplot(312)
plot([Ps(:) Ps(:)])
subplot(313)
plot([Pcs(:) Pcs(:)])







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% I could estiamte with general routine (Nov. 18, 2011)
% =>  Estimate alpha, beta, Q, and R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% [alpha beta Vrev Q R] = deal(params(1,end), params(2,end), params(3,end), params(4,end), params(5,end));

%% initialize
alpha = alpha0;
beta = beta0;
Vrev = Vrev0;
Q = Q0;
R = R0;
Ke = Ke0(:);

% fix Ke
Ke = Ke0;

alpha = rand();
beta = rand();

mses = [NaN];
LLs = [NaN];


paramTrue = [alpha0; beta0; Vrev0; Q0; R0];
params = [alpha; beta; Vrev;  Q; R];


for itr = 1:20
    %% 1. Kalman smoothe with wrong param (params)

%     Ke = Kes(:,end)'
    [Xs Ps Pcs] = kalman_smth_1d(Y-Vrev, Ke*II, alpha, beta, 1, 1, Q, R, Xo, Po);
    mses(itr+1) = (X-Xs)*(X-Xs)'/ns 
    
        
    %% 2. M-step: update parameter and calculate LL

    [alpha, beta, ~, ~, Q, R, LL] = m_step_abcd(Y-Vrev0, Ke0*II, Xs, reshape(Ps,1,1,[]), reshape(Pcs,1,1,[]), alpha, beta, 1, 1, Q, R)
    
    
    params(1,itr+1) = alpha;
    params(2,itr+1) = beta;
    params(3,itr+1) = Vrev;
    params(4,itr+1) = Q;
    params(5,itr+1) = R;
    
    
    LLs(itr+1) = LL
    
    
    
    
    
    clf;
    subplot(411)
    plot(Y); title ('measurement')
    box off
    subplot(412)
    plot([X' Xs'])
    title (sprintf('V_m (mse=%.3f)',mses(itr+1)))
    box off
    subplot(425)
    plot([paramTrue(1:5,end) params(1:5,end)], '.-'); title ('param')
    box off
    subplot(426)
    plot(params'); title ('param over iteration')
    box off
    subplot(427)
    plot(LLs)
    box off
    subplot(428)
    plot(mses)
    title ('mse')
    box off

%     subplot(427)
%     plot([KeTrue Kes(:,end)],'.-')
    drawnow    
end

saveas(1,'EstimLinGeneral_fixed_Ke.pdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Next, estimate Ke as well. (Nov. 18, 2011)
% There is redundancy in parameterization that general routine cannot distinguish
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = alpha0;
B = beta0*Ke0;
C = 1;
D = Ke0;
Q = Q0;
R = R0;

Ke = Ke0;

paramTrue = [A Q R]';

alpha = 0.1

params=[A Q R];
Bs=[B(:)];
Ds=[D(:)];
mses = [NaN];
LLs = [NaN];
for itr = 1:20
    %% 1. Kalman smoothe with wrong param (params)
    
    [Xs Ps Pcs] = kalman_smth_1d(Y-Vrev, II, A, B, 1, D, Q, R, Xo, Po);
    mses(itr+1) = (X-Xs)*(X-Xs)'/ns 
    
        
    %% 2. M-step: update parameter and calculate LL

    [A, B, ~, D, Q, R, LL] = m_step_abcd(Y-Vrev0, II, Xs, reshape(Ps,1,1,[]), reshape(Pcs,1,1,[]), A, B, 1, D, Q, R)
    
    
    params(1,itr+1) = A;
    params(2,itr+1) = Q;
    params(3,itr+1) = R;
    
    
    Bs(:,itr+1) = B(:);
    Ds(:,itr+1) = D(:);
    LLs(itr+1) = LL
    
    
    
    
    %%
    numRow = 5;
    numCol = 2;
    clf;
    subplot(numRow,numCol,1)
    plot(Y); title ('measurement')
    box off
    subplot(numRow,numCol,2)
    plot([X' Xs'])
    title (sprintf('V_m (mse=%.3f)',mses(itr+1)))
    box off
    subplot(numRow,numCol,3)
    plot([paramTrue(1:3,end) params(1:3,end)], '.-'); title ('param')
    box off
    subplot(numRow,numCol,4)
    plot(params'); title ('param over iteration')
    box off
    subplot(numRow,numCol,5)
    plot([beta0*Ke0(:)  B(:)]); title ('estimated B')
    box off
    subplot(numRow,numCol,6)
    plot(Bs'); title ('estimated B over iteration')
    box off
    
    subplot(numRow,numCol,7)
    plot([ Ke0(:) D(:)]); title ('estimated D')
    box off
    subplot(numRow,2,8)
    plot(Ds'); title ('estimated B over iteration')
    box off
    
    subplot(numRow,numCol,9)
    plot(LLs); title ('LL')
    box off
    subplot(numRow,numCol,10)
    plot(mses); title ('mse')
    box off

%     subplot(427)
%     plot([KeTrue Kes(:,end)],'.-')
    drawnow    
end
    

saveas(1,'EstimLinGeneral_unconstrained.pdf')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare with Cpp version!

clear
clc
load Y.txt -ascii;
load II.txt -ascii;
load Xs.txt -ascii;
load Ps.txt -ascii;
load Pcs.txt -ascii;
load params.txt -ascii;
load Ke.txt -ascii;



disp ('initial param')
[alpha, beta, Vrev, Q, R]  = deal(params(1),params(2),params(3),params(4),params(5))

% M-step by Matlab 
disp ('Call M-step by Matlab')
[Mat.alpha, Mat.beta, Mat.Vrev, Mat.Q, Mat.R, Mat.Ke, Mat.LL] = m_step_full_1d_new_param(Y, II, Xs, Ps, Pcs, alpha, beta, Vrev, Q, R, Ke);

% M-step by C
disp ('Call M-step by Cpp')
[cpp.alpha, cpp.beta, cpp.Vrev, cpp.Q, cpp.R, cpp.K, cpp.LL] = m_step_1d(Y, II(1,:), Xs, Ps, Pcs, alpha, beta, Vrev, Q, R, Ke);

Mat
cpp