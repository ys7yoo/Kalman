% separated into a function 
% originally based on test3_LDS_paramEstim_keFix.m
%
% Test fitting of model parameters but with Ke fixed.
%
% X = latent state
% I = injected current (filtered and added to X)
% Y = observed state 
% U = ke * I  % filtered injected current
%
% Basic equations
% X_t+1 = a*X_{t} + beta*U_t + noiseX  
% Y_t = X_t + mu_y + U_t + noiseY
%
% X ~ N(D*U*beta, Q)
% Y|X ~ N(X+muy+U, R),   R = sigy^2*I
% Y ~ (N(D*U*beta+U+muy, Q + sig^2*I)
%
% Note: from linear regime, we plan to 
% KEEP: ke, R, beta?
% RE-ESTIM: a, Q, mu_y?, 


function [prshat] = estimParamGradientDescent(yy, U, ntaps, prs0)



% addpath jptools
% addpath Kalman/jptools



% %% Compute marginal likelihood by hand (only tractable for small # samples)
% % Y ~ (N(D*U*beta+U+muy, Q + sig^2*I)
% Qy = inv(Qinv)+sigy^2*I;
% resid = (yy-xmu-muy-U);
% normresid = norm(resid);
% trm1 = -.5*logdet(Qy);
% trm2 = -.5*resid'*inv(Qy)*resid;
% ll = trm1+trm2



%%  4. Optimize likelihood for LDS parameters
% ---------------------------------
% Y|prs ~ N(U, C)
% 
%  cov C = inv(Q) + sigy^2*I
% mean U = beta*Miapp*ke + U + muy
% --------------------------------------


%% Perform optimization (fminunc)
lFunc = @(prs)(neglogmargLi_LDSvec_noKe(prs,yy,ntaps,U));

tic;
opts = optimset('display', 'iter','largescale','off','maxiter',2000,'maxfuneval',25000);
prshat = fminunc(lFunc,prs0, opts);
toc;

% %% Perform optimization again (fmincon);  % UNCOMMENT TO RUN FMINCON VERSION
% prs1 = prs0; nprs = length(prs1);
% prs1(ntaps+(1:2)) = exp(prs1(ntaps+(1:2)));
% LB = zeros(nprs,1);
% LB(1:ntaps) = -ntaps; UB(1:ntaps) = ntaps;
% LB(ntaps+(1:2))=0.001;UB(ntaps+(1:2))=1000;
% 
% LB(ntaps+3)=-100;
% UB = inf(nprs,1);
% lFunc2 = @(prs)(neglogmargLi_LDSvec2(prs,yy,ntaps,II));
% 
% tic;
% opts = optimset('display', 'iter','largescale','off','algorithm','active-set','maxiter',200,'maxfuneval',2500);
% prshat2 = fmincon(lFunc2,prs1,[],[],[],[],LB,UB,[],opts);
% toc;





return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Routine for generating data (Originally in the beginning of
%% test3_LDS_paramEstim_keFix.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
% 1. Set up LDS simulation 
% ns =  5000; % number of samples
ns = input ('number of samples? ')
sigIinj = 1; % stdev of injected current
sigxmarg = .2; % marginal stdev of noise in latent state
sigy = .1;    % observation noise
muy = -50;   % mean of y 
beta = .01;   % multiplier on current's effect on Vm (should be << 1)

% Model parameters (params to be estimated later)
rr = [.9 .3 .2]';  % specify poles (roots) of IIR filter
avec = real(poly(rr))';  % coefficients of IIR filter
ntaps = length(avec)-1;

% Find the prior marginal covariance of the AR process
%  (solution to: C = A*C*A' + Q )
acov = compARautocorr(-avec(2:end),1);
sigx = sigxmarg/sqrt(acov(1)); % noise added in each time step
acov = acov*sigx.^2; % autocovariance

C = toeplitz(acov);  % prior covariance matrix
% Check for positive-definiteness
if min(eig(C))<0
    warning('non-stationary AR process! C has negative eigenvalues');
end

% Make Electrode filter
ke = [.07 .1 .15 .25 .65 1.4  .25 .05]';
nke = length(ke);
plot(ke);


%% 2.  Generate sample from LDS model
%  -----------------------------------

% Make injected current
Iapp = randn(ns,1)*sigIinj; % applied current
II = [0;Iapp(1:end-1)]; % current gets injected on *next* timestep
U = filter(flipud(ke),1,II); % compute current filtered by Ke

xmu = beta*filter(1,avec,U); % Noiseless membrane response
%xmu2 = filter(flipud(ke),avec,beta*II); % Noiseless membrane response (equiv)

% Make noisy membrane response
xplusns = beta*U + randn(ns,1)*sigx;
xx = filter(1,avec,xplusns);  % latent variable (filter with dynamics filter)
yy = xx + muy + U + randn(ns,1)*sigy;  % Sample from the observed variable

% Make plot
subplot(221);
plot([beta*U,xplusns]); title('noisy current injected');
subplot(222);
plot(beta*U,xplusns,'.'); title('raw vs. noisy current');
axis equal; xlabel('current'); ylabel('current + noise');
subplot(223);
plot([xx+muy,yy]); title('latent variable and measurement');
subplot(224);
plot(xx+muy,yy,'.'); title('obs vs. latent variable');
axis equal; xlabel('x'); ylabel('y');


%% 3. Run Kalman Smoother (with true parameters): estimate X|Y
%  ------------------------------------------------------------

% Compute mean 
xmu = beta*filter(flipud(ke),avec,II); % mean

% Compute inverse covariance
acovvec = compARautocorr(-avec(2:end),sigx^2,ntaps*2+1); % autocovariance
Qi = inv(toeplitz(acovvec));  % inverse covariance (small piece)
Qd = spdiags(Qi,-ntaps:ntaps);
Qdd = [Qd(1:ntaps,:); repmat(Qd(ntaps+1,:),ns-2*ntaps,1);Qd(ntaps+2:end,:)];
Qinv = spdiags(Qdd,-ntaps:ntaps,ns,ns);

% MAP estimate for X|Y,Iinj  (KS estimate)
I = speye(ns);
xhat = (Qinv+I/sigy^2)\((yy-muy-U)/sigy^2 + Qinv*xmu); % Kalman smoothing in one line!

% --------------------------
% make plot
subplot(211);
t = (1:ns)';
h = plot(t,xmu,'b', t,yy-muy-U,'g',t,xhat,'r',t,xx,'k');
%set(h(3:4), 'linewidth', 1.5);
legend('Xprior', 'Y', 'Xmap', 'Xtrue','location','northeast');

subplot(223);
xl = min(xx):.1:max(xx);
plot(xl,xl,'k',xx,xmu,'bo', xx,yy-muy-U,'go', xx,xhat, 'r.'); axis tight;
xlabel('true X'); ylabel('estim X');

subplot(224);
plot(xl, xl*0,'k', xx,xx-xmu, 'bo', xx,xx-yy+muy+U,'go', xx,xx-xhat,'r.'); axis tight;
xlabel('true X'); ylabel('error');

fprintf('\nAvg Errors:\n----------\n');
fprintf(' Xprior: %.2f\n      Y: %.2f\n   Xmap: %.2f\n', [std(xx-xmu) std(xx-yy) std(xx-xhat)]);
% ---------------------------





%% Pick initial values
muy0 = mean(yy); % estimate of mean added to y
ke0 = .2*ones(nke,1)+.1;
sigx0 = 2;
sigy0 = 2;
avec0 = zeros(ntaps,1);

% Initial parameters (transformed nonlinearly to make fitting easier)
prstrue = [-avec(2:end); log(1./sigx.^2); log(1./(sigy^2)); muy; beta];
prs0 =  [avec0; log(1./sigx0.^2); log(1./(sigy0^2)); muy0; beta*2];
%prs0 = prstrue*.75;

U = filter(flipud(ke),1,II);
Ltrue = -neglogmargLi_LDSvec_noKe(prstrue,yy,ntaps,U);% log-likelihood at true params 
L0 = -neglogmargLi_LDSvec_noKe(prs0,yy,ntaps,U); % log-likelihood at initial params

fprintf('logli of True params (%.2f), Initial params (%.2f)\n',Ltrue,L0);







%%  4. Optimize likelihood for LDS parameters
% ---------------------------------
% Y|prs ~ N(U, C)
% 
%  cov C = inv(Q) + sigy^2*I
% mean U = beta*Miapp*ke + U + muy
% --------------------------------------

% call this function 
prshat = estimParamGradientDescent(yy, U, ntaps, prs0)



%% --- Report true vs. estimated parameters ----
ahat = prshat(1:ntaps);
sigxhat = sqrt(exp(-prshat(ntaps+1)));
sigyhat = sqrt(exp(-prshat(ntaps+2)));
muyhat = prshat(ntaps+3);
betahat = prshat(ntaps+4:end);
fprintf('\nParam comparison: true / estim\n');
fprintf('true avec:  ');
for j = 1:ntaps
    fprintf('%.2f, ', -avec(1+j));
end
fprintf('\nestim avec: ');
for j = 1:ntaps
    fprintf('%.2f, ', ahat(j));
end
fprintf('\n');
fprintf('sigx:   %.3f  /   %.3f \n', sigx, sigxhat);
fprintf('sigy:   %.3f  /   %.3f \n', sigy, sigyhat);
fprintf('muy:  %.3f  / %.3f \n', muy, muyhat);
fprintf('beta: %.3f  / %.3f \n', beta, betahat);

% % ----- Make plot -------
% subplot(211);
% tk = 1:nke;
% plot(1,beta,'ko', 1,,'k--',tk,betahat, 'r-o'); axis tight;
% legend('true k', 'initial k', 'ML k estim', 'location', 'northeast');
% title('ke');

subplot(211);
tprs = 1:ntaps;
plot(tprs,prstrue(tprs),tprs,prs0(tprs),'k--',tprs,prshat(tprs),'r');
title('recursive filter (taps)');

subplot(212);
tprs = [ntaps+1:ntaps+2, ntaps+4];
plot(tprs,prstrue(tprs),tprs,prs0(tprs),'k--',tprs,prshat(tprs),'ro-');
legend('true prs', 'initial prs', 'ML prs', 'location', 'northeast');
title('other params (sigx, xigy, beta)');








%% save results
%save multiTapJP
save multiTapJp.mat avec beta sigx sigy xx yy U muy sigxmarg

%% run multiple times 

%ns = input('number of samples? ');

ns = 50000
%for rep = 1:10
for rep = 11:100
    test3vec_LDS_paramEstim_keFix
    filename = sprintf('multiTapJp_N%d_rep%d.mat',ns,rep)
    save(filename, 'avec', 'beta', 'sigx', 'sigy', 'xx', 'yy', 'U', 'muy', 'sigxmarg','ntaps','prshat')
end


