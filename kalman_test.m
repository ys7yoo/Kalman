



%% set parameters 
clear
dt = .01;
T = 200;
Rep = 100;


% paremeters for neuron & measurement model
gL = .1;
EL = -65;
A = 1-gL*dt;
B = dt;
C = 1;
D = 1;


%A = [1-gL*dt gL*dt; 0 1];        
%B = [dt; 0];
%B = [0; 0];
%C = [1 0];  
%D = [0; 0];

% generate data for multiple level of noise
Q = 1;
R = 1;

Xo = EL;
Po = 1;
    
% control input 
Idc = 0;
Isig = 10;
Itau = 5/dt;
%for itr = 1:Rep        
I = Isig*randn(1,T) + Idc;

% Generate data
[X,Y] = generate_lds(I, A, B, C, D, Q, R, Xo, Po);

%
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
    
SAVE_PARAM_AND_DATA = 0
if (SAVE_PARAM_AND_DATA)
    %% save 1 dim version
    save ('A1.txt', 'A', '-ascii')
    save ('B1.txt', 'B', '-ascii')
    save ('C1.txt', 'C', '-ascii')
    save ('D1.txt', 'D', '-ascii')
    save ('Q1.txt', 'Q', '-ascii')
    save ('R1.txt', 'R', '-ascii')
    save ('Xo1.txt', 'Xo', '-ascii')
    save ('Po1.txt', 'Po', '-ascii')
    
    save ('X1.txt', 'X', '-ascii');
    save ('Y1.txt', 'Y', '-ascii');
    save ('U1.txt', 'I', '-ascii');
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% filtering v.s. smoothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call smoother. Then smoother calls filter internally.
[xsmooth Ps Pcs xfilt Pf Pcf] = kalman_smth(Y, I, A, B, C, D, Q, R, Xo, Po);
sprintf('mse_filter: %.2f',mean(squeeze(Pf(1,1,:))))
sprintf('mse_smoother: %.2f',mean(squeeze(Ps(1,1,:))))

%
V = X(1,:)';
Vf = xfilt(1,:)';
Vs = xsmooth(1,:)';

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
Bnd = Q*5;
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


%% change 1-dim vesion for KF
% ori.xfilt =xfilt;
% ori.Pf = Pf;

% run new version
disp('run KF ori')
clear gen;
tic;
[gen.Xp gen.Pp gen.Xf gen.Pf gen.Kn gen.LL] = kalman_filt(Y, I, A, B, C, D, Q, R, Xo, Po);
toc;

disp('run KF 1D')
clear one;
tic;
[one.Xp one.Pp one.Xf one.Pf one.Kn one.LL] = kalman_filt_1d(Y, I, A, B, C, D, Q, R, Xo, Po);
toc;

max(abs(one.Xp - gen.Xp))
max(abs(one.Pp' - squeeze(gen.Pp)))
max(abs(one.Xf - gen.Xf))
max(abs(one.Pf' - squeeze(gen.Pf)))
max(abs(one.Kn - gen.Kn))
max(abs(one.LL - gen.LL))

% => more than two times fast. The problem was how to access memory!



%%
disp('run KS general')
clear gen;
tic;
[gen.xsmooth gen.Ps gen.Pcs gen.xfilt gen.Pf gen.LL] = kalman_smth(Y, I, A, B, C, D, Q, R, Xo, Po);
toc;

disp('run KS 1d')
clear one;
tic;
%[one.xsmooth one.Ps one.Pcs one.xfilt one.Pf one.LL] = kalman_smth_1d(Y, I, A, B, C, D, Q, R, Xo, Po);
[one.xsmooth one.Ps one.Pcs] = kalman_smth_1d(Y, I, A, B, C, D, Q, R, Xo, Po);  % test mex version
toc;

max(abs(one.xsmooth - gen.xsmooth))
max(abs(one.Ps' - squeeze(gen.Ps)))
max(abs(one.Pcs(2:end)' - squeeze(gen.Pcs(2:end))))

clf
subplot(311)
plot([one.xsmooth'   gen.xsmooth'])
subplot(312)
plot([one.Ps'   squeeze(gen.Ps)])
subplot(313)
plot([one.Pcs(2:end)'   squeeze(gen.Pcs(2:end))])
% max(abs(one.xfilt - gen.xfilt))
% max(abs(one.Pf' - squeeze(gen.Pf)))
% max(abs(one.LL - gen.LL))

% => OK. two times faster!
% run KS general
% Elapsed time is 0.029667 seconds.
% run KS 1d
% Elapsed time is 0.011824 seconds.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare with general version kalman_smth_gen.mex (from kalman_smth_gen.cpp)



disp('run KS general')
clear gen;
tic;
[gen.xsmooth gen.Ps gen.Pcs gen.xfilt gen.Pf gen.LL] = kalman_smth(Y, I, A, B, C, D, Q, R, Xo, Po);
toc;

disp('run KS 1d')
clear one;
tic;
%[one.xsmooth one.Ps one.Pcs one.xfilt one.Pf one.LL] = kalman_smth_1d(Y, I, A, B, C, D, Q, R, Xo, Po);
[gen.xsmooth gen.Ps gen.Pcs] = kalman_smth_gen(Y, I, A, B, C, D, Q, R, Xo, Po);  % test mex version
toc;

max(abs(one.xsmooth - gen.xsmooth))
max(abs(one.Ps' - squeeze(gen.Ps)))
max(abs(one.Pcs(2:end)' - squeeze(gen.Pcs(2:end))))

clf
subplot(311)
plot([one.xsmooth'   gen.xsmooth'])
subplot(312)
plot([one.Ps'   squeeze(gen.Ps)])
subplot(313)
plot([one.Pcs(2:end)'   squeeze(gen.Pcs(2:end))])





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now estimate parameters using emKalman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc

% test 1. start from true value
paramInit.A = A;
paramInit.B = B;
paramInit.C = C;
paramInit.D = D;
paramInit.Q = Q;
paramInit.R = R;
paramInit.Xo = Xo;
paramInit.Po = Po;
[paramEstim] = em_kalman_abcd(Y, I, paramInit); % general version
%[paramEstim LLs] = em_kalman(Y, UU, paramInit, EM, paramTrue) % electrod compensation version

paramInit

paramEstim


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's take a closer look at each step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_iter = 50;






