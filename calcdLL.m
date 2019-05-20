function [dLL] = calcdLL(Xs, Y, U, Ps, Pcs, Theta)

N = length(Y);                         
M=size(U,1);
[alpha beta gamma Ke Q R Xo Po] = deal(Theta(1), Theta(2), Theta(3), ...
                                       Theta(4:M+3), Theta(M+4), Theta(M+5), ...
                                       Theta(M+6), Theta(M+7)); 
Ke = Ke';
                                   
                                   
%% now let's calc derivatives 
C=1;
S1 = sum(Xs(2:end));
S0 = S1+Xs(1)-Xs(end);
SU0 = sum(U(:,1:end-1)')';

S = Xs*Xs' + sum(Ps,3);
S11 = S - (Xs(:,1)*Xs(:,1)'+ Ps(:,:,1));
S10 = Xs(:,2:end)*Xs(:,1:end-1)' + sum(Pcs(:,:,2:N),3);
S00 = Xs(:,1:end-1)*Xs(:,1:end-1)' + sum(Ps(:,:,1:N-1),3);

SXU = Xs*U';
SXU00 = SXU - Xs(:,end)*U(:,end)';
SXU10 = Xs(:,2:end)*U(:,1:end-1)';
SUU = U*U';
SUU00 = SUU - U(:,end)*U(:,end)';

SYU = Y*U';
SYY = Y*Y'+C*sum(Ps,3)*C';
SYX = Y*Xs';


%% calculate derivatives 
%Ke = (SYU-C*SXU)/SUU;      % simple approximation
Qi = 1/Q;   % inv(Q);
Ri = 1/Q;   % inv(R);
dLdABC = Qi*([S10 SXU10*Ke' S1]-[alpha beta gamma]*[S00 SXU00*Ke' S0; Ke*SXU00' Ke*SUU00*Ke' Ke*SU0; S0 SU0'*Ke' N-1]);
dLdKe = (beta*Qi*(SXU10-alpha*SXU00-gamma*SU0')+Ri*(SYU-C*SXU)) - Ke*(Qi*beta.^2*SUU00+Ri*SUU);

ABCnew = [S10 SXU10*Ke' S1]/[S00 SXU00*Ke' S0; Ke*SXU00' Ke*SUU00*Ke' Ke*SU0; S0 SU0'*Ke' N-1];
dLdQ = 0.5*((S11-ABCnew*[S10'; Ke*SXU10'; S1])*Qi-(N-1))*Qi;

dLdR = 0.5*((SYY + C*S*C' + Ke*SUU*Ke'- Ke*SYU' - SYU*Ke' - C*SYX' -SYX*C'  ...
             + C*SXU*Ke' + Ke*SXU'*C')*Ri -N)*Ri;


%% update initial condition
dLdXo = (Xs(:,1)-Xo)/Ps(:,:,1);
dLdPo = 0.5*((Xs(:,1)-Xo)*(Xs(:,1)-Xo)'/Po-1)/Po;

dLL = [dLdABC dLdKe dLdQ dLdR dLdXo dLdPo];
%dLL = [dLdABC; dLdKe; dLdQ; dLdR; dLdXo; dLdPo];
