function [alpha beta gamma Ke Q R Xo Po]=m_step_full(alpha, beta, gamma, Ke, Q, R, Y, U, Xo, Po, Xs, Ps, Pcs, constraint)


%% FIXME: Xo and Po are included in Xs and Ps
if (nargin<14)
    constraint = 0;
end

% M-step for EMKalman
% simplified for 1D input

N = size(Y,2);
dimX = size(alpha,1);  % => should be 1

%% M-step: update parameters

C=1;
S1 = sum(Xs(2:end));
S0 = S1+Xs(1)-Xs(end);
SU0 = sum(U(:,1:end-1),2);

SS = Xs*Xs';
%SS11 = SS - Xs(:,1)*Xs(:,1)';
SS10 = Xs(:,2:end)*Xs(:,1:end-1)';
SS00 = SS - Xs(:,end)*Xs(:,end)';


KeU = Ke*U;

P = sum(Ps,3);
P11 = P - Ps(:,:,1);
P10 = sum(Pcs(:,:,2:N),3);
P00 = P - Ps(:,:,end);

% % SP =  SS + P;
% % SP11 = SS11 + P11;
% % SP10 = SS10 + P10;
% % SP00 = SS00 + P00;

SXU = Xs*U';
SXU00 = SXU - Xs(:,end)*U(:,end)';
SXU10 = Xs(:,2:end)*U(:,1:end-1)';
SUU = U*U';
SUU00 = SUU - U(:,end)*U(:,end)';

SYU = Y*U';
% % SYY = Y*Y'; %+C*sum(Ps,3)*C';
% % SYX = Y*Xs';

%% calculate new Ke
Qi = inv(Q);
Ri = inv(R);
%Ke = (SYU-C*SXU)/SUU;      % simple approximation
KeNew = (beta*Qi*(SXU10-alpha*SXU00-gamma*SU0')+Ri*(SYU-C*SXU))/(Qi*beta.^2*SUU00+Ri*SUU);

%% calculate new alpha, beta, gamma
%     % update alpha and beta
%     AB = [SP10 SXU10*Ke']/[SP00 SXU00*Ke'; Ke*SXU00' Ke*SUU00*Ke'];      %AB = SP10*inv(SP00);
%     alpha = AB(:,1:dimX);
%     beta = AB(:,dimX+1:end);
%ABC = [SP10 SXU10*Ke' S1]/[SP00 SXU00*Ke' S0; Ke*SXU00' Ke*SUU00*Ke' Ke*SU0; S0 SU0'*Ke' N-1];
% Different implementation (Jul. 26, 2010)
XKeU0 = Xs(:,1:end-1)*KeU(1:end-1)';
SKeU0 = sum(KeU(1:end-1));
ABC = [SS10+P10 Xs(:,2:end)*KeU(1:end-1)' S1] / [SS00+P00 XKeU0 S0; XKeU0 KeU(1:end-1)*KeU(1:end-1)' SKeU0; S0 SKeU0' N-1];





%% update Q
%Q=(SP11-AB*[SP10'; Ke*SXU10'])/(N-1);   % equivalently Q=(SP11-SP10*inv(SP00)*SP10')/(N-1);
%Q=(SP11-ABC*[SP10'; Ke*SXU10'; S1])/(N-1);   % equivalently Q=(SP11-SP10*inv(SP00)*SP10')/(N-1);
% Different implementation (Jul. 26, 2010)
temp=(Xs(2:end)-alpha*Xs(1:end-1)-beta*KeU(1:end-1)-gamma); % same form as L1
Q = (temp*temp' + P11 -2*alpha*P10 + alpha^2*P00)/(N-1);    
Q=(Q+Q')/2;


%% update R
% % R = (SYY + C*SP*C' + Ke*SUU*Ke'- Ke*SYU' - SYU*Ke' - C*SYX' -SYX*C'  ...
% %      + C*SXU*Ke' + Ke*SXU'*C')/N;
% Different implementation (Jul. 26, 2010)
R=((Y-Xs-KeU)*(Y-Xs-KeU)' + P )/N;  % same form as L2
R=diag(diag(R));        % diagonalize
% % 
% % if (R<=0)
% %     LL = -NaN;
% %     paramEstim = [];
% %     return;
% % end


%% update alpha, beta, gamma, and Ke here
alpha = ABC(:,1:dimX);
beta = ABC(:,dimX+1:dimX+1);
gamma = ABC(:,dimX+2:end);

%     if alpha > 1
%         alpha = .99;
%     end

if (constraint)
    %Ke = (KeNew>1e-3).*KeNew;    % non-negative constraint
    Ke = max(KeNew, 0);   
else
    Ke = KeNew;
end


%% update initial condition
Xo=Xs(:,1);
Po=Ps(:,:,1);



% % %% calc LL with update parameters 
% % KeU = Ke*U;
% % 
% % temp=(Xs(2:end)-alpha*Xs(1:end-1)-beta*KeU(1:end-1)-gamma);
% % L1 = temp*temp' + P11 -2*alpha*P10 + alpha^2*P00;
% % L1 = L1 / Q;
% % L2 = (N-1)*log(abs(Q));
% % 
% % L3 = (Y-Xs-KeU)*(Y-Xs-KeU)' + P;
% % L3 = L3/R;
% % L4 = N*log(abs(R));
% % 
% % L5 = (Xs(1)-Xo).^2/Po;
% % L6 = log(abs(Po));
% % 
% % LL = -0.5*(L1+L2+L3+L4+L5+L6);


