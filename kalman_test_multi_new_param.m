%% set parameters 
clear

dt = .01;
T = 5000;
Rep = 100;


% paremeters for neuron & measurement model
gL = .1;
EL = -65;
%A = [0.7 0.2 0.1; 1 0 0; 0 1 0]
A = [0.7 0.1 0.05; 1 0 0; 0 1 0]
B = [0.1; 0; 0]
C = [ 1 0 0]
D = 1


%A = [1-gL*dt gL*dt; 0 1];        
%B = [dt; 0];
%B = [0; 0];
%C = [1 0];  
%D = [0; 0];

% generate data for multiple level of noise
Q = [0.1 0 0; 0 0 0; 0 0 0];
R = 0.1;

Xo = zeros(3,1);
Po = eye(3,3);
    
% control input 
Idc = 0;
Isig = 10;
Itau = 5/dt;
%for itr = 1:Rep        
I = Isig*randn(1,T) + Idc;

%% Generate data
[X,Y] = generate_lds(I, A, B, C, D, Q, R, Xo, Po);
% X = X + EL;
% Y = Y + EL;

%% Or, load from the file 





%% 
clf;
subplot(311)
plot (dt*(1:T),I)
xlabel('ms')
ylabel('\muA')
title ('Input current')
subplot(312)
plot (dt*(1:T),X(1,:))
xlabel('ms')
ylabel('mV')
title (sprintf('Membrane voltage (Q=%.0e)',Q(1,1)))
%     subplot(413)
%     plot (dt*(1:T),Po(:,:,idxParam))
%     xlabel('ms')
%     ylabel('mV')
%     title (sprintf('Reverse potential (Q1=%f)',Q0))
subplot(313)
plot (dt*(1:T),Y)
xlabel('ms')
ylabel('mV')
title (sprintf('Measured voltage (R=%.0e)',R))

drawnow
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% filtering v.s. smoothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call smoother. Then smoother calls filter internally.
[xsmooth Ps Pcs LL xfilt Pf Pcf] = kalman_smth(Y, I, A, B, C, D, Q, R, Xo, Po);
sprintf('mse_filter: %.2f',mean(squeeze(Pf(1,1,:))))
sprintf('mse_smoother: %.2f',mean(squeeze(Ps(1,1,:))))

%
V = X(1,:)';
Vf = reshape(xfilt(1,1,:),[],1);
Vs = reshape(xsmooth(1,1,:),[],1);

clf                                                   
subplot (3,1,1)
plot(V, 'r+:'); hold on
plot(Vf, 'ms--');
plot(Vs, 'bo-');
legend ('V', 'X_{filt}', 'X_{smooth}')
title (sprintf('State estimate (Q=%1.1e,R=%.0f)',Q,R))
subplot (3,1,2)
errFilt=V-Vf;
errSmth=V-Vs;
plot (errFilt, 'ms-');hold on
plot (errSmth, 'bo-'); 
plot (squeeze(Pf(1,1,:)), 'rs:');hold on
plot (squeeze(Ps(1,1,:)), 'cs:');hold on

%legend (sprintf('filter (mse = %.2f)', var(x'-xfilt')),sprintf('smoother (mse = %.2f)', var(x'-xsmooth')))
legend ('filter', 'smoother', 'Pf', 'Ps')
title ('Error')
subplot (3,3,7)
%B= ceil(max(abs([errFilt;errSmth])));
Bnd = max(Q(:)*5);
bin = [-Bnd:Bnd/10:Bnd];
hist (errFilt,bin); hold on
set(gca,'xlim', [-Bnd Bnd]);
title (sprintf('Histogram of filtering error \n(mean=%.2f, var=%.2f)', mean(errFilt), var(errFilt)))
subplot (3,3,8)
hist (errSmth,bin); hold on
set(gca,'xlim', [-Bnd Bnd]);
title (sprintf('Histogram of smoothing error \n(mean=%.2f, var=%.2f)', mean(errSmth), var(errSmth)))
subplot (3,3,9)
plot (X(1,:)',Vf, 'ms'); hold on
plot (X(1,:)',Vs, 'bo');
m = min([xfilt(1,:),xsmooth(1,:)]);
M = max([xfilt(1,:),xsmooth(1,:)]);
plot ([m M], [m M], 'k--')


xlabel ('X_{true}')
ylabel ('Estimate')
legend('filter', 'smoother')

axis equal



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compare filtering of Cpp with matlab (Feb. 1, 2012)

dim = size(A,1)
% load data
load YY.txt
load UU.txt 

% Kf result by Cpp
load Xf.txt
XfCpp = Xf;
load Pf.txt
PfCpp = Pf;
%PfCpp = reshape(Pf,dim,dim,[]);
load Xp.txt
XpCpp = Xp;
load Pp.txt
PpCpp = Pp;
%PpCpp = reshape(Pp,dim,dim,[]);

% Kf by Matlab
[Xp Pp Xf Pf Kf LL] = kalman_filt(YY, UU, A, B, C, D, Q, R, Xo, Po)

Pp = reshape(Pp,dim*dim,[]);
Pf = reshape(Pf,dim*dim,[]);

%% compare state estimate
clf
subplot(321)
plot(XpCpp')
subplot(323)
plot(Xp')
subplot(325)
plot(XpCpp'-Xp')


subplot(322)
plot(XfCpp')
subplot(324)
plot(Xf')
subplot(326)
plot(XfCpp'-Xf')
% =>  numerical error

%% compare variance estimate
clf
idx=1:10; % good
idx=length(Y)-10:length(Y);

subplot(321)
plot(PpCpp(1,idx)')
subplot(323)
plot(Pp(1,idx)')
subplot(325)
plot(PpCpp(1,idx)'-Pp(1,idx)')


subplot(322)
plot(PfCpp(1,idx)')
subplot(324)
plot(Pf(1,idx)')
subplot(326)
plot(PfCpp(1,idx)'-Pf(1,idx)')
% =>  Check 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compare Kalman smoother of Cpp with that of matlab (Feb. 1, 2012)

dim = size(A,1)
% load data
load YY.txt
load UU.txt 

% Kf result by Cpp
load Xs.txt
XsCpp = Xs;
load Ps.txt
PsCpp = Ps;
%PfCpp = reshape(Pf,dim,dim,[]);
load Pcs.txt
PcsCpp = Pcs;
%PpCpp = reshape(Pp,dim,dim,[]);

% Kf by Matlab
tic
[Xs Ps Pcs] = kalman_smth(YY, UU, A, B, C, D, Q, R, Xo, Po);
toc

Ps = reshape(Ps,dim*dim,[]);
Pcs = reshape(Pcs,dim*dim,[]);

%% compare state/variance estimate
%idx=length(Y)-10:length(Y);
%idx=1:10;
idx = 1:length(Y);
clf
subplot(331)
plot(XsCpp(:,idx)')
subplot(334)
plot(Xs(:,idx)')
subplot(337)
plot(XsCpp(:,idx)'-Xs(:,idx)')


subplot(332)
plot(PsCpp(:,idx)')
subplot(335)
plot(Ps(:,idx)')
subplot(338)
plot(PsCpp(:,idx)'-Ps(:,idx)')

subplot(333)
plot(PcsCpp(:,idx)')
subplot(336)
plot(Pcs(:,idx)')
subplot(339)
plot(PcsCpp(:,idx)'-Pcs(:,idx)')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now estimate parameters using emKalman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
A
B


% test 1. start from true value
paramInit.A = A;% + randn(size(A));
paramInit.B = B;% + randn(size(B));
paramInit.C = C;
paramInit.D = D;
paramInit.Q = Q;
paramInit.R = R;
paramInit.Xo = Xo;
paramInit.Po = Po;
paramInit.Vrev = 0;

EM.max_iter = 50;
EM.checkConverged = 0;
EM.checkDecreased = 1;
EM.eps = 1e-5;
    
updateParam = [1 1 0 0 1 1 0];  % which parameter to update? A B C D Q R Vrev
[paramEstim] = em_kalman_abcd_new_param(Y, I, paramInit, EM, updateParam); % general version

% 
% paramInit
% 
% paramEstim

disp ('estimated params (A, B, C, D, Q, R)')
paramEstim.A
paramEstim.B
% paramEstim.C
% paramEstim.D
paramEstim.Q
paramEstim.R


%% plot result
clf;
subplot(321)
plot(paramEstim.LLs)


subplot(323)
plot(paramEstim.A(1,:),'o-'); hold on
plot(A(1,:),'r+--')
title ('the 1st row of A')

subplot(324)
eigAest = eig(paramEstim.A)
plot(real(eigAest),imag(eigAest),'o'); hold on
eigA = eig(A)
plot(real(eigA),imag(eigA),'r+'); hold on
axis equal
title ('eigen values of A')
legend ('estimate','true')
box off


pHat = [paramEstim.B(1,1); paramEstim.Q(1,1); paramEstim.R(1,1)]
pTrue = [B(1,1); Q(1,1); R(1,1)]

subplot(325)
plot(pHat,'o'); hold on;
plot(pTrue,'+r')
legend ('estimate','true')
set(gca,'xtick',1:3)
set(gca,'xticklabel',{'B','Q','R'})

box off

% old interface 
% Ainit = rand(size(A));
% Binit = rand(size(B));
% Dinit = rand(size(D));
%Dinit = D;
% [Aest, Best, Cest, Dest, Qest, Rest, Xoest, Poest, LL] = ...
%     em_kalman(Y, I, Ainit, Binit, C, Dinit, Q, R, Xo, Po, max_iter);
% Aest
% Best
% Dest
% Qest
% Rest
% clf
% semilogy(LL)
% 

%% compare estimated state with true
clf;
subplot(311)
plot(Y)
title ('measurement')

subplot(312)
plot(paramEstim.Xs(1,:)); hold on
plot(V,'r--')
legend('estimate','true')
title (sprintf('estimated state (MSE=%.2f)',(paramEstim.Xs(1,:)-V')*(paramEstim.Xs(1,:)'-V)/length(V)))



return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compile KS to speed up
ns = length(Y);
emlmex kalman_smth.m -eg {emlcoder.egs(0,[1,ns]) emlcoder.egs(0,[1,ns]) emlcoder.egs(0,size(A)) emlcoder.egs(0,size(B)) emlcoder.egs(0,size(C)) emlcoder.egs(0,size(D))  emlcoder.egs(0,size(Q)) emlcoder.egs(0,size(R)) emlcoder.egs(0,size(Xo)) emlcoder.egs(0,size(Po))}




