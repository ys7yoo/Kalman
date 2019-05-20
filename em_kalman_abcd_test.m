


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% test more general version (Oct. 16, 2011)
% % clear
% % dimX = 3;
% % 
% % A = randn(dimX,dimX)
% % eig(A)
% % B = randn(dimX,1)
% % Q = 0.1*eye(dimX,dimX)
% % 
% % % one-dim measurement => Not so good
% % C = randn(1,dimX)
% % D = 0
% % R = 0.1
% % 
% % % measurement dimension is increased 
% % C = randn(dimX,dimX)
% % D = randn(dimX,1)
% % R = 0.1*eye(dimX,dimX)




% save params

%% simpler param with dim = 3
clear
clc

dimX = 3;
A = randn(dimX,dimX)
[U D V] = svd(A);
% replace singular value
D = diag(rand(dimX,1))
A = U*D*V';
eig(A)


B = randn(dimX,1)
Q = 0.01*eye(dimX,dimX)

dimY = 3;
dimU = 1;
C = randn(dimY,dimX)
D = randn(dimY,dimU)
R = 0.01*eye(dimY,dimY)


% %% simpler param with dim = 1
% clear
% clc
% 
% dimX = 1;
% A = 0.9;
% eig(A)
% 
% 
% B = randn(dimX,1)
% Q = 0.01*eye(dimX,dimX)
% 
% C = randn(dimX,dimX)
% D = randn(dimX,1)
% R = 0.01*eye(dimX,dimX)


%% generate data using LDS
Xo = zeros(dimX,1)
Po = 0.01*eye(dimX,dimX)


ns = input('number of samples? ');
if isempty(ns)
    ns = 5000;
end

U = randn(1,ns);
U = filter([0.25 0.25 0.25 0.25], [1 0.9], U);
U = filter([0.25 0.25 0.25 0.25], 1, U);
U = filter([0.25 0.25 0.25 0.25], 1, U);
U = filter([0.25 0.25 0.25 0.25], 1, U);
U = filter([0.25 0.25 0.25 0.25], 1, U);
%U = sin((1:ns)/25)+ 
%U = ones(1,ns);


[X,Y] = generate_lds(U, A, B, C, D, Q, R, Xo, Po);


xrange = [0 1000];
clf
subplot(311)
plot(U')
title ('input')
set(gca,'xlim', xrange)

subplot(312)
plot(X')
title ('state')
set(gca,'xlim', xrange)
%
subplot(313)
plot(Y')
title ('measurement')
set(gca,'xlim', xrange)

% save data
 


%% compile KS to speed up
PRE_COMPILE = 0
if (PRE_COMPILE)
    %%
    emlmex kalman_smth.m -eg {emlcoder.egs(0,[dimY,ns]) emlcoder.egs(0,[dimU,ns]) emlcoder.egs(0,size(A)) emlcoder.egs(0,size(B)) emlcoder.egs(0,size(C)) emlcoder.egs(0,size(D))  emlcoder.egs(0,size(Q)) emlcoder.egs(0,size(R)) emlcoder.egs(0,size(Xo)) emlcoder.egs(0,size(Po))}
end

%% Kalman smoothing
disp ('Kalman smooting')
tic
[Xs Ps Pcs LL] = kalman_smth(Y, U, A, B, C, D, Q, R, Xo, Po);
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
plot(Xs')
title ('estimated state')
set(gca,'xlim', xrange)
%
subplot(515)
plot(Xs'-X')
title ('error')
set(gca,'xlim', xrange)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now, estimate parameters 


%EM.num_trial = 1;
EM.num_trial = 8;
EM.max_iter = 20;
EM.checkConverged = 1;
EM.checkDecreased = 1;
EM.eps = 1e-6;


clf; 
clear paramInit paramEst paramEstim LLtrace;


% initialize 
for trial = 1:EM.num_trial
    paramInit(trial).A = A + 0.1*rand(size(A));
    paramInit(trial).B = B + 0.1*rand(size(B)); 
    paramInit(trial).C = C ;
    paramInit(trial).D = D + 0.1*rand(size(D));
    paramInit(trial).Q = Q + 0.01*diag(rand(dimX,1)); 
    paramInit(trial).R = R + 0.01*diag(rand(dimY,1)); 
    paramInit(trial).Xo = Xo;   %X(1) + sqrt(5)*randn();
    paramInit(trial).Po = Po;   % 5*rand();
    
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
    [paramEstim(trial)] = em_kalman_abcd(Y, U, paramInit(trial), EM, [1 1  0 1 1 1]);
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
LLend  = [paramEstim(:).LL];
[LLmax idx] = max(LLend);
paramEst = paramEstim(idx);


%% plot estimatie parameters 
clf
row=3;
col=4;
subplot(row,col,1); hold on;
for trial = 1:length(paramEstim)
    plot(paramEstim(trial).LLs);
end
plot(paramEst.LLs,'r'); % plot the best one

subplot(row,col,5)
plot (A(:), 'xr'); hold on
plot (paramEst.A(:),'o')
%set(gca, 'ylim', [0 1.1])
title ('A')

%plot eig values
subplot(row,col,9)
eigA = eig(A)
plot(real(eigA),imag(eigA),'xr'); hold on
eigAest = eig(paramEst.A)
plot(real(eigAest),imag(eigAest),'o'); hold on
axis equal
title ('eigen values of A')
%legend ('estimate','true')
box off



subplot(row,col,6)
plot (B(:), 'xr-'); hold on
plot (paramEst.B(:),'o-')
title ('B')

subplot(row,col,7)
plot (C(:), 'xr-'); hold on
plot (paramEst.C(:),'o-')
title('C')

subplot(row,col,8)
plot (D(:), 'xr-'); hold on
plot (paramEst.D(:),'o-')
title ('D')

subplot(row,col,11)
plot (diag(Q), 'xr-'); hold on
plot (diag(paramEst.Q),'o-')
title ('diag(Q)')

subplot(row,col,12)
plot (diag(R), 'xr-'); hold on
plot (diag(paramEst.R),'o-')
title ('diag(R)')


return

saveas(1,'fitting_multi_dim.pdf')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EM kalman with membrane kernel Ke
%% Kernel in both the state equation and measurement equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% run like T=10000;max_iter=100;em_kalman_abc_test
tic

%T=20000;
T=1000;
max_iter=100;

%% generate data

addpath ../YsTools
%clear
%T = 1000;
dt = .1;

% electrode kernel
Ke = .5*[0 0 .8 7 3 .3 0];
M = length(Ke);

% parameters
A = 1-.1*.1;
B = [.1/sum(Ke)*Ke -.65];
C = 1;
D = [Ke 0];

Q = 1e0;
R = 1e0;

% initial conditions
Xo = -80;
Po = 1;

% control input 
Idc = 0;
Isig = 10;
Itau = 5/dt;
%for itr = 1:Rep        
Iapp = Isig*randn(1,T) + Idc;
% store history terms to each column
II=stackRows(Iapp,M,1);

%% generate data
IIext = [II; ones(1,T)];


[X,Y] = generate_lds(IIext, A, B, C, D, Q, R, Xo, Po);


%
clf;
subplot(311)
plot (dt*(1:T), Iapp);
%plot (dt*(1:T),[Iinj(:) Ifil(:)])
%legend ('I_{inj}', 'I_{filterd}')
xlabel('ms')
ylabel('nA')
title ('Input current')
set(gca,'xlim',[0 100])
subplot(312)
plot (Ke,'x-')
subplot(313)
plot (dt*(1:T),Y(:),'k'); hold on;
plot (dt*(1:T),X(1,:),'b--'); hold on
plot (dt*(1:T),(Ke*II),'r--');
legend ('Vrec', 'Vm', 'Ve')
xlabel('ms')
ylabel('mV')
title ('Membrane and electrode voltage')
set(gca,'xlim',[0 100])
%title (sprintf('Membrane and electrode voltage (Q=%.0e)',Q(1,1)))
% subplot(413)
% plot (dt*(1:T),Y(:))
% xlabel('ms')
% ylabel('mV')
% title (sprintf('Measured voltage (Vm+Ve+noise)(R=%.0e)',R))

filename = sprintf('generated_Q%.0f_R%.0f_T%d_itr%d.png', Q,R,T,max_iter);
saveas(1, filename)

drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check the performance of smoother
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Xs Ps Pcs LL] = kalman_smth(Y, IIext, A, B, C, D, Q, R, Xo, Po);

%
clf
plot (X', '+-r'); hold on
plot (Xs', 'o-')
legend ('true', 'estimated')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate A,beta,D,Q,R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EM.num_trial = 1;
%EM.num_trial = 8;
EM.max_iter = 100;

% Initial guess of kernel using simple linear regression
Klin = Y'\II';
Klin = flipud(Klin);
Klin = Klin.*5.55./sum(Klin);

clf; 
clear paramInit paramEst paramEstim LLtrace;
LLend = [];
XXs = zeros(EM.num_trial,T);
parfor trial = 1:EM.num_trial
    
    %% trial 
    disp(sprintf('%dth trial',trial))
    
    % initialize
    paramInit(trial).A = 0.9 + 0.1*rand();
    paramInit(trial).B = [Klin + randn(1,M)  -rand()];
    paramInit(trial).C = C;
    paramInit(trial).D = [Klin/sum(Klin) + randn(1,M) 0];%rand(size(Ke));
    paramInit(trial).Q = Q*(.5+rand(size(Q)));
    paramInit(trial).R = R*(.5+rand(size(R)));
    paramInit(trial).Xo = X(1) + sqrt(5)*randn();
    paramInit(trial).Po = 5*rand();
    
%     % start some of the value from the original value (for debug)
%     paramInit(trial).A = A;
%     paramInit(trial).B = B;
%     paramInit(trial).C = C;
%     paramInit(trial).D = D;
%     paramInit(trial).Q = Q;
%     paramInit(trial).R = R;
%     paramInit(trial).Xo = Xo;
%     paramInit(trial).Po = Po;
       
    % estimate parameters
    [paramEstim(trial)] = em_kalman_abcd(Y, IIext, paramInit(trial), EM);
    LLend = [LLend paramEstim(trial).LL];
    
    disp(sprintf('  LL=%.2e', LL(end)));
    subplot(121); hold on;
    semilogy(paramEstim(trial).LLs);
    %title ('Log likelihood')
    title ('L(\Theta)')
    xlabel ('Iteration')
    subplot(122); hold on;
    plot(trial, LL,'x');
    title ('Final Log likelihood')
    xlabel ('Trial')
    drawnow;
end
saveas(1, 'EM_Kalman_ABCD_trials.pdf')

% choose parameter set with maximum LL
[LLmax idx] = max(LLend);
paramEst = paramEstim(idx);


%% check accuracy of estimatie parameters 
clf
subplot(221)
plot (A, 'xr'); hold on
plot (paramEst.A,'o')
set(gca, 'ylim', [0 1.1])
title ('A')
subplot(222)
plot (B(1:end-1), 'xr-'); hold on
plot (paramEst.B(1:end-1),'o-')
title ('B(1:end-1)')
subplot(223)
plot (B(end), 'xr-'); hold on
plot (paramEst.B(end),'o-')
set(gca, 'ylim', [1.2*min(B(end),paramEst.B(end)) 0])
title('B(end)')
subplot(224)
plot (D, 'xr-'); hold on
plot (paramEst.D,'o-')
title ('D')

saveas(1, 'EM_Kalman_ABCD_params.pdf')

return;

%%
Ifil = filter(Ke, 1, Iinj);
IfilEst = filter(paramEst.Ke, 1, Iinj);

% convert back to original parameters
C = dt/beta/sum(Ke);
gL = (1-A)/beta/sum(Ke);
Vo = gamma ./ (1-A);

paramEst.C = dt/paramEst.beta/sum(paramEst.Ke);
paramEst.gL = (1-paramEst.A)/paramEst.beta/sum(paramEst.Ke);
paramEst.Vo = paramEst.gamma ./ (1-paramEst.A);


% save result
filename = sprintf('estimation_Q%.0f_R%.0f_T%d_itr%d.mat', Q,R,T,max_iter);
save(filename)


%% draw figures
subplot(612);hold off
plot (dt*(1:T), Iinj,'k--'); hold on;
plot (dt*(1:T), Ifil,'r');
plot (dt*(1:T), IfilEst);
%plot (dt*(1:T),[Iinj(:) Ifil(:)])
legend ('I_{app}', 'True I_{inj}', 'Estimated I_{inj}')
xlabel('ms')
ylabel('\muA')
title ('Input current')
subplot(613);hold off
plot (dt*(1:T),paramEst.Xs(1,:),'r'); hold on
plot (dt*(1:T),X(1,:))
legend ('True V_m', 'Estimated V_m')
xlabel('ms')
ylabel('mV')
title (sprintf('Membrane voltage (Q=%.0e)',Q(1,1)))
subplot(614)
plot (dt*(1:T),Ke*II,'r'); hold on
plot (dt*(1:T),paramEst.Ke*II)
legend ('True V_e', 'Estimated V_e')
xlabel('ms')
ylabel('mV')
title (sprintf('Electrode voltage (Q=%.0e)',Q(1,1)))
% compare A
subplot(6,4,17);hold off
plot (A, '+r'); hold on
plot (paramEst.A, 'o')
%legend ('True', 'Estimated')
title ('\A')
set(gca,'ylim', [min([0; paramEst.A]) 2*max([A; paramEst.A])]);
% compare beta
subplot(6,4,18);hold off
plot ([beta], '+r'); hold on
plot ([paramEst.beta], 'o')
%legend ('True', 'Estimated')
title ('\beta')
set(gca,'ylim', [min([0; paramEst.beta]), 2*max([beta; paramEst.beta])]);
% compare gamma
subplot(6,4,19);hold off
plot ([gamma], '+r'); hold on
plot ([paramEst.gamma], 'o')
%legend ('True', 'Estimated')
title ('\gamma')
set(gca,'ylim', [2*min([gamma; paramEst.gamma]), max([0; paramEst.gamma]), ]);
% compare kernel
subplot(6,4,20); hold off
plot (Ke', '+-r'); hold on
%plot (Klin, 'x-');
plot (paramEst.Ke, 'o-');
%legend ('True', 'Estimated')
title ('Ke')
%
% original parameter
subplot(6,4,21);hold off
plot (C, '+r'); hold on
plot (paramEst.C, 'o')
%legend ('True', 'Estimated')
title ('C')
set(gca,'ylim', [min([0; paramEst.C]), 2*max([C; paramEst.C])]);
% compare gL
subplot(6,4,22);hold off
plot (gL, '+r'); hold on
plot (paramEst.gL, 'o')
%legend ('True', 'Estimated')
title ('g_L')
set(gca,'ylim', [min([0; paramEst.gL]) 2*max([gL; paramEst.gL])]);
% compare gamma
subplot(6,4,23);hold off
plot (Vo, '+r'); hold on
plot (paramEst.Vo, 'o')
%legend ('True', 'Estimated')
title ('Vo')
set(gca,'ylim', [2*min([Vo; paramEst.Vo]), max([0; paramEst.Vo]), ]);
% compare kernel
subplot(6,4,24);hold off
plot (Ke', '+-r'); hold on
plot (paramEst.Ke, 'o-');
legend ('True', 'Estimated')
title ('Ke')


%%
filename = sprintf('estimation_Q%.0f_R%.0f_T%d_itr%d.png', Q,R,T,max_iter);
saveas(1, filename)

%%
toc

return

%% draw separate figures

%% 0) input and output
clf;
subplot(211)
plot (dt*(1:T), Iinj);
%plot (dt*(1:T),[Iinj(:) Ifil(:)])
%legend ('I_{inj}', 'I_{filterd}')
xlabel('ms')
ylabel('nA')
title ('Input current')
set(gca,'xlim',[0 100])
subplot(212)
plot (dt*(1:T),Y(:),'k'); hold on;
plot (dt*(1:T),X(1,:),'b--'); hold on
plot (dt*(1:T),(Ke*II),'r--');
legend ('Vrec', 'Vm', 'Ve')
xlabel('ms')
ylabel('mV')
title ('Measured voltage')
set(gca,'xlim',[0 100])
%title (sprintf('Membrane and electrode voltage (Q=%.0e)',Q(1,1)))
% subplot(413)
% plot (dt*(1:T),Y(:))
% xlabel('ms')
% ylabel('mV')
% title (sprintf('Measured voltage (Vm+Ve+noise)(R=%.0e)',R))

filename = sprintf('generated_Q%.0f_R%.0f_T%d_itr%d.png', Q,R,T,max_iter);
saveas(1, filename)


%% 1) draw LL traces
clf
hold on;
for trial = 1:length(LLtrace)
    semilogy(real(LLtrace(trial).LLs));
end
title ('L(\Theta)')
xlabel ('Iteration')

%% 2) Parameter estimation of LDS
disp('Estimated parameters for LDS')
disp(sprintf('A=%.3f', paramEst.A))
disp(sprintf('beta=%.3f', paramEst.beta))
disp(sprintf('gamma=%.3f', paramEst.gamma))
% compare A
clf
subplot(141);hold off
plot (A, '+r'); hold on
plot (paramEst.A, 'o')
%legend ('True', 'Estimated')
title ('\A')
set(gca,'ylim', [min([0; paramEst.A]) 2*max([A; paramEst.A])]);
% compare beta
subplot(142);hold off
plot ([beta], '+r'); hold on
plot ([paramEst.beta], 'o')
%legend ('True', 'Estimated')
title ('\beta')
set(gca,'ylim', [min([0; paramEst.beta]), 2*max([beta; paramEst.beta])]);
% compare gamma
subplot(143);hold off
plot ([gamma], '+r'); hold on
plot ([paramEst.gamma], 'o')
%legend ('True', 'Estimated')
title ('\gamma')
set(gca,'ylim', [2*min([gamma; paramEst.gamma]), max([0; paramEst.gamma]), ]);
% compare kernel
subplot(144); hold off
plot (Ke', '+-r'); hold on
%plot (Klin, 'x-');
plot (paramEst.Ke, 'o-');
legend ('True', 'Estimated')
title ('Ke')

%% 3) Kernel only

paramEst.A
clf
plot (Ke', '+-r', 'MarkerSize',10); hold on
%plot (Klin, 'x-');
plot (paramEst.Ke, 'o-',  'MarkerSize',10);
legend ('True', 'Estimated')
title ('Ke')
set(gca,'ylim', [-0.1 4])

%% 4) Original parameters 
disp('Estimated parameters for LDS')
disp(sprintf('C=%.3f', paramEst.C))
disp(sprintf('g_L=%.3f', paramEst.gL))
disp(sprintf('V_o=%.3f', paramEst.Vo))

% original parameter
clf
subplot(141);hold off
plot (C, '+r'); hold on
plot (paramEst.C, 'o')
%legend ('True', 'Estimated')
title ('C')
set(gca,'ylim', [min([0; paramEst.C]), 2*max([C; paramEst.C])]);
% compare gL
subplot(142);hold off
plot (gL, '+r'); hold on
plot (paramEst.gL, 'o')
%legend ('True', 'Estimated')
title ('g_L')
set(gca,'ylim', [min([0; paramEst.gL]) 2*max([gL; paramEst.gL])]);
% compare gamma
subplot(143);hold off
plot (Vo, '+r'); hold on
plot (paramEst.Vo, 'o')
%legend ('True', 'Estimated')
title ('Vo')
set(gca,'ylim', [2*min([Vo; paramEst.Vo]), max([0; paramEst.Vo]), ]);
% compare kernel
subplot(144);hold off
plot (Ke', '+-r'); hold on
plot (paramEst.Ke, 'o-');
legend ('True', 'Estimated')
title ('Ke')

%% Estimation results
clf
disp('mse of Ve')
Ve = Ke*II;
VeEst = paramEst.Ke*II;
mse1=(Ve-VeEst)*(Ve-VeEst)'/length(Ve)
disp('mse of Vm')
mse2=(X-Xs)*(X-Xs)'/length(Xs)


subplot(211)
plot (dt*(1:T),Ke*II,'r'); hold on
plot (dt*(1:T),paramEst.Ke*II)
legend ('True V_e', 'Estimated V_e')
xlabel('ms')
ylabel('mV')
title (sprintf('Electrode voltage (MSE=%.4f)',mse1))
set(gca,'xlim',[0 100])
subplot(212);hold off
plot (dt*(1:T),X(1,:));hold on
plot (dt*(1:T),paramEst.Xs(1,:),'r');
legend ('True V_m', 'Estimated V_m')
xlabel('ms')
ylabel('mV')
title (sprintf('Membrane voltage (MSE=%.2f)',mse2))
set(gca,'xlim',[0 100])

%% draw figures
subplot(612);hold off
plot (dt*(1:T), Iinj,'k--'); hold on;
plot (dt*(1:T), Ifil,'r');
plot (dt*(1:T), IfilEst);
%plot (dt*(1:T),[Iinj(:) Ifil(:)])
legend ('I_{app}', 'True I_{inj}', 'Estimated I_{inj}')
xlabel('ms')
ylabel('\muA')
title ('Input current')


%% to compare with previous results

EM = paramEst
save est X Y II A beta C Ke Q R Xo Po EM




