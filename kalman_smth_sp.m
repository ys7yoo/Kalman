function [xhat postVar2 xmu varOff] = kalman_smth_sp(u, yy, param)

%% 2. Run Kalman Smoother (with true parameters): estimate X|Y
%  ------------------------------------------------------------

% Compute prior mean of X | Iinj
xmu = filter(1,[1;-param.alpha],[param.xo; u(1:end-1)]);     % prediction without dynamics noise


% Form Qinv (inverse covariance matrix of X).
ns = length(yy);
I = speye(ns);
Qdiag = spdiags([0;ones(ns-2,1);0],0,ns,ns); % identity with corners removed
Qoff = spdiags(ones(ns,2),[-1 1], ns,ns); % +1/-1 diagonals
Qinv = (I+param.alpha.^2*Qdiag - param.alpha*Qoff)/param.sigx.^2; % inverse covariance matrix

% MAP estimate for X|Y,Iinj
xhat = (Qinv+I/param.sigy^2)\((yy-param.muy-param.beta*u)/param.sigy^2 + Qinv*xmu); % Kalman smoothing in one line!


if nargout >1
    Linv = Qinv+I/param.sigy^2; % inverse posterior covariance
%     tic;
%     postVar1 = diag(inv(Linv)); % the explicit inverse of the posterior cov
%     toc
%     tic
    postVar2 = invtridiag_special(Linv(1,1),Linv(2,1),Linv(2,2),ns); % fast version
%    toc

end

if nargout >3
    % variance of offdiagonal 
    varOff = diag(inv(Linv),1);
end
    


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test routine 

clc
tic
[xhat postVar xmu varOff] = kalman_smth_sp(u, yy, paramLDS);
toc

expectedErrorFromPosterio = full(mean(postVar))

var(xhat-xx)
    
%
% --------------------------
% make plot
clf
subplot(211);
t = (1:ns)';
h = plot(t,xmu,'b', t,yy-paramLDS.muy,'g',t,xhat,'r',t,xx,'k');
%set(h(3:4), 'linewidth', 1.5);
legend('Xprior', 'Y', 'Xmap', 'Xtrue','location','northeast');

subplot(223);
xl = min(xx):.1:max(xx);
plot(xl,xl,'k',xx,xmu,'bo', xx,yy-paramLDS.muy,'go', xx,xhat, 'r.'); axis tight;
xlabel('true X'); ylabel('estim X');

subplot(224);
plot(xl, xl*0,'k', xx,xx-xmu, 'bo', xx,xx-yy+paramLDS.muy,'go', xx,xx-xhat,'r.'); axis tight;
xlabel('true X'); ylabel('error');

fprintf('\nAvg Errors:\n----------\n');
fprintf(' Xprior: %.2f\n      Y: %.2f\n   Xmap: %.2f\n', [std(xx-xmu) std(xx-yy) std(xx-xhat)]);
% ---------------------------


%% Compare with YS's KS 
addpath('../Kalman')

Y = yy'-paramLDS.muy;
%U = Miapp';
U = u';     % filtered input (Ke*Iapp)
A = paramLDS.alpha;
%B = ke';
B = 1;
C = 1;
%D = zeros(1,length(ke));
D = paramLDS.beta;;
Q = paramLDS.sigx^2;
R = paramLDS.sigy^2;
Xo = paramLDS.xo;
Po = paramLDS.sigxmarg^2;
tic
[Xs Ps Pcs LL Xf Pf Xp Pp] = kalman_smth_1d(Y, U, A, B, C, D, Q, R, Xo, Po);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kalman smoothing for the folling model
% X(k+1) = A*X(k) + B*U(k) + v(k)
% Y(k) = C*X(k) + D*U(k) + w(k)
% Ev = 0, Evv' = Q
% Ew = 0, Eww' = R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% compare result
clf
t = (1:ns)';
subplot(411)
plot(t,xx,'+k-'); hold on
%errorbar (t, xhat, sqrt(postVar2),'bo')
%errorbar (t, Xs, sqrt(Ps),'rs')
plot(t, xhat, 'ro--')
plot(t, Xs, 'bs:')
title ('estimated mean')
legend ('true', 'JP', 'YS');


subplot(412)
%plot(t,xx,'+k-'); hold on
%errorbar (t, xhat, sqrt(postVar2),'bo')
%errorbar (t, Xs, sqrt(Ps),'rs')
plot(t, xhat-xx, 'ro--'); hold on
plot(t, Xs-xx', 'bs:')
title (sprintf('error (MSE.JP=%.2f,MSE.YS=%.2f)',[var(xhat-xx) var(Xs-xx')]))
legend ('JP', 'YS');


subplot(413)
plot(t, postVar, 'ro--'); hold on
plot(t, Ps, 'bs:')
legend ('JP', 'YS');
title (sprintf('estimated var (JP=%.4f,YS=%.4f)',[full(mean(postVar)) mean(Ps)]))


subplot(414)
plot(t(2:end), varOff, 'ro--'); hold on
plot(t(2:end), Pcs(2:end), 'bs:')
legend ('JP', 'YS');
title (sprintf('estimated lag-one var (E[x_k x_{k-1}] (JP=%.4f,YS=%.4f)',[full(mean(varOff)) mean(Pcs(2:end))]))


%% Make sure they give the exactly same result!
clf
% JP result
subplot(321);
xl = min(xx):.1:max(xx);
plot(xl,xl,'k', xx,xhat, 'r.'); axis tight; axis equal
xlabel('true X'); ylabel('estim X');
title ('estimate from JP')


subplot(322);
plot(xl, xl*0,'k', xx,xhat-xx,'r.'); axis tight; axis equal
xlabel('true X'); ylabel('error');
title (sprintf('error from JP (mean=%.2f,std=%.2f)',[mean(xhat-xx) std(xhat-xx)]))

% YS result
subplot(323);
xl = min(xx):.1:max(xx);
plot(xl,xl,'k', xx,Xs, 'b.'); axis tight; axis equal
xlabel('true X'); ylabel('estim X');
title ('estimate from YS')

subplot(324);
plot(xl, xl*0,'k', xx,Xs'-xx,'b.'); axis tight; axis equal
xlabel('true X'); ylabel('error');
title (sprintf('error from YS (mean=%.2f,std=%.2f)',[mean(Xs'-xx) std(Xs'-xx)]))


subplot(325)
plot (xl,xl,'k',xhat, Xs, 'o')
xlabel('JP')
ylabel('YS')
title ('mean')
axis equal

subplot(326)
pl = min(full(min(postVar)),min(Ps)):0.05:max(full(max(postVar)),max(Ps));
plot (pl,pl, 'k',postVar, squeeze(Ps), 'o')
xlabel('JP')
ylabel('YS')
title ('variance')
axis equal



fprintf('\nMSE of JP:\n----------\n');
fprintf(' Xprior: %.2f\n      Y: %.2f\n   Xmap: %.2f\n expected:%.2f\n', [var(xx-xmu) var(xx-yy) var(xx-xhat) expectedErrorFromPosterio]);
fprintf('\nMSE of YS:\n----------\n');
fprintf(' smth: %.2f\n expected: %.2f\n', [var(xx-Xs') mean(Ps)]);
