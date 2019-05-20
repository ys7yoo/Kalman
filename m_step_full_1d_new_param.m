function [alpha beta Vrev Q R Ke LL] =m_step_full_1d_new_param(Y, U, Xs, Ps, Pcs, alpha, beta, Vrev, Q, R, Ke, EM)
% update parameters given measurement (Y), input (U), and Kalman smoother results (Xs, Ps, Pcs) starting from initial param (paramInit)


N = size(Y,2);
% nke = size(KeInit,1);


% make sure Ke is a row vector
Ke = reshape(Ke,1, []);

if (nargin<12)
     constraint = 0;
     repeat = 1;
     fixBeta = 0;
else 
	constraint = EM.MstepConstraint;
    repeat = EM.MstepRepeat;
    fixBeta = EM.fixBeta;
end






%% M-step: update parameters
C=1;

% calc square sums
SS00 = Xs(1:end-1)*Xs(1:end-1)';    % sum Xn*Xn' from n=1 to N-1
SS  = SS00 + Xs(end)*Xs(end)';      % sum Xn*Xn from n=1 to N
SS11 = SS - Xs(:,1)*Xs(:,1)';       % sum Xn*Xn from n=2 to N  
SS10 = Xs(:,2:end)*Xs(:,1:end-1)';  % sum X{n+1}*xn' from n=1 to N-1

% calc terms from variance 
P = sum(Ps);
P11 = P - Ps(1);
P10 = sum(Pcs(2:N));
P00 = P - Ps(end);

% terms involves with the input
SXU = Xs*U';
SXU00 = SXU - Xs(end)*U(:,end)';
SXU10 = Xs(:,2:end)*U(:,1:end-1)';
SUU = U*U';
SUU00 = SUU - U(:,end)*U(:,end)';

% terms involves with the measurement
Ysub = Y-Vrev;
SYY = Ysub*Ysub';
SYX = Ysub*Xs';
SYU = Ysub*U';
%% I can calculate SYU using only the 1st row of U
I=U(1,:);
nke = length(Ke);
(Y-Vrev)*I'
(Y(2:end)-Vrev)*I(1:end-1)'
% the last term
(Y(nke:end)-Vrev)*I(1:end-nke+1)'

%% calculate new Ke

KeU = Ke*U;
% this is also equal to  
% filter(Ke, 1, U(1,:));

Qi = 1/Q;
Ri = 1/R;

%Ke = (SYU-C*SXU)/SUU;      % simple approximation for very small beta << 1
KeNew = (beta*Qi*(SXU10-alpha*SXU00)+Ri*(SYU-C*SXU))/(Qi*beta.^2*SUU00+Ri*SUU);

%% calculate new alpha, beta with the new parameterization (Aug. 17, 2011)
XKeU10 = SXU10*Ke';
XKeU00 = SXU00*Ke';
AB = [SS10+P10 XKeU10] / [SS00+P00 XKeU00; XKeU00 Ke*SUU00*Ke'];



%% update Q
quad1 = (SS11+P11) + alpha*(SS00+P00)*alpha + beta*Ke*SUU00*Ke'*beta ...
       -2*alpha*(SS10+P10) + 2*alpha*SXU00*Ke'*beta - 2*SXU10*Ke'*beta;
Q = quad1/(N-1);    
%     Q=(Q+Q')/2;


%% update R

quad2 = SYY + C*(SS+P)*C' + Ke*SUU*Ke'- Ke*SYU' - SYU*Ke' - C*SYX' -SYX*C'  ...
             + C*SXU*Ke' + Ke*SXU'*C';        
R = quad2/N;
%     R=diag(diag(R));        % diagonalize


%% update Vrev
Vrev = mean(Y-Xs-KeU);



%% update alpha, beta, Vrev, and Ke here
alpha = AB(1);      %AB(:,1:dimX);

if ~fixBeta
    beta = AB(2);       %AB(:,dimX+1:dimX+1);
end

%     if alpha > 1
%         alpha = .99;
%     end
Ke = KeNew;



%% calculate the log likelihood
L1 = trace(quad1)/Q + (N-1)*log(Q);
L2 = trace(quad2)/R + N*log(R);
LL = -0.5*(L1+L2);





