function L = neglogmargLi_LDSvec_noKe(prs,y,ntaps,II)
% L = neglogmargLi_LDSvec_noKe(prs,y,ntaps,II)
%
% Negative log-marginal likelihood of parameters for a vector linear-dynamical
% system (i.e., Kalman filtering model).  For use in ML fitting.
%
% INPUT:
%  prs = parameters:
%        avec - AR(p) coefficients
%        rhox - log of inverse noise variance (precision) in dynamics in each time step
%        rhoy - log of inverse noise variance (precision) in observations
%        muy - mean added to observations y
%        beta - log of scalar on the effect of injected current on X (latent)
%  y = observations (vector)
%  II = injected current (zero-shifted by 1 bin).
%
% OUTPUT: 
%  L = negative marginal log-likelihood
%
% --- Model equations: -----
% X_t+1 = a*X_{t} + beta*U_t + noiseX  
% Y_t = X_t + mu_y + U_t + noiseY
%
% X ~ N(D*U*beta, Q)                  % latent variable
% Y|X ~ N(X+mu+U, R), R = sigy^2*I    % conditional observed
% Y ~ (N(D*U*beta+U+mu, Q + R)  % marginal observed


% 1. Unpack parameters
avec = [1; -prs(1:ntaps)];
rhox = exp(prs(ntaps+1)); % marginal variance of noise in X (latent)
rhoy = exp(prs(ntaps+2)); % variance of noise in Y|X 
muy = prs(ntaps+3);
beta = exp(prs(ntaps+4)); % multiplier on linput to latent variable
slen = length(II);

% Compute inverse covariance estimate
acovvec = compARautocorr(-avec(2:end),1/rhox,ntaps*2+1); % autocovariance

% Check stability of LDS
C = toeplitz(acovvec(1:ntaps));
if min(eig(C))<1e-15
    %warning('Cannot evaluate likelihood: unstable roots of LDS equation');
    L = 1e25; % arbitrary large value
else
% ----- compute log-likelihood ----------

% 1. Filter injected current with ke & dynamics params
U = filter(1,avec,II); % mean
ydiff = y-beta*U-muy-II; % mean-centered Y

% 2. Form inverse covariance matrices for X|Iinj and Y|X
Qi = inv(toeplitz(acovvec));  % inverse covariance (small piece)
Qd = spdiags(Qi,-ntaps:ntaps);
Qdd = [Qd(1:ntaps,:); repmat(Qd(ntaps+1,:),slen-2*ntaps,1);Qd(ntaps+2:end,:)];
Qx = spdiags(Qdd,-ntaps:ntaps,slen,slen);  % inverse cov for X
Qy = Qx/rhoy+speye(slen); % marginal covariance for Y divided by Qx^-1 (convenient representation)

% 4. Compute the marginal log-likelihood
logDetTrm = .5*(logdet(Qy)-logdet(Qx));  % log determinant term
LambdaInvY = (Qy\(Qx*ydiff));  % needed for quadratic term
QuadTrm = .5*ydiff'*LambdaInvY;  % quadratic term
L = logDetTrm+QuadTrm;

end