function [LL, dLL] = calcNegLLdLL(Xs, Y, U, Ps, Pcs, Xo, Po, Theta)

% The goal is to calculate negative log likelihood, -log P(Y|Theta)
% given Y, Theta, and P(X|Y,Theta) (Xs, Ps, Pcs), 



N = length(Y);                         
M=size(U,1);

[alpha beta gamma Ke Q R] = deal(Theta(1), Theta(2), Theta(3), ...
                                       Theta(4:M+3), Theta(M+4), Theta(M+5));

Ke = Ke';



%% reimplementation (Jun. 25, 2010)

% calc basic terms.
    
KeU = Ke*U;
% general dimension version
% P = sum(Ps,3);
% P11 = P - Ps(:,:,1);
% P10 = sum(Pcs(:,:,2:N),3);
% P00 = P - Ps(:,:,end);
% one dimensional version
P = sum(Ps);
P11 = P - Ps(1);
P10 = sum(Pcs(2:N));
P00 = P - Ps(end);

%% calc LL
temp=(Xs(2:end)-alpha*Xs(1:end-1)-beta*KeU(1:end-1)-gamma);
L11 = temp*temp' + P11 -2*alpha*P10 + alpha^2*P00;
L1 = L11 / Q;
L2 = (N-1)*log(abs(Q));

L33 = (Y-Xs-KeU)*(Y-Xs-KeU)' + P;
L3 = L33/R;
L4 = N*log(abs(R));

L5 = (Xs(1)-Xo).^2/Po;
L6 = log(abs(Po));

LL = .5*(L1+L2+L3+L4+L5+L6);


if nargout >1 % calc gradient too
    Qi = 1/Q;
    Ri = 1/R;
    
    C = 1;
    S1 = sum(Xs(2:end));
    S0 = S1+Xs(1)-Xs(end);
    SU0 = sum(U(:,1:end-1),2);

    SS = Xs*Xs';
    %SS11 = SS - Xs(:,1)*Xs(:,1)';
    SS10 = Xs(:,2:end)*Xs(:,1:end-1)';
    SS00 = SS - Xs(:,end)*Xs(:,end)';
    

    SXU = Xs*U';
    SXU00 = SXU - Xs(:,end)*U(:,end)';
    SXU10 = Xs(:,2:end)*U(:,1:end-1)';
    SUU = U*U';
    SUU00 = SUU - U(:,end)*U(:,end)';

    SYU = Y*U';
    % % SYY = Y*Y'; %+C*sum(Ps,3)*C';
    % % SYX = Y*Xs';

    %% dL/dABC


    XKeU0 = Xs(:,1:end-1)*KeU(1:end-1)';
    SKeU0 = sum(KeU(1:end-1));
    A1 = [SS10+P10 Xs(:,2:end)*KeU(1:end-1)' S1];
    A2 = [SS00+P00 XKeU0 S0; XKeU0 KeU(1:end-1)*KeU(1:end-1)' SKeU0; S0 SKeU0' N-1];
    dLdABC = 2*Qi*( A1 - [alpha beta gamma]*A2);
    
    
    %% calculate new Ke
    %Ke = (SYU-C*SXU)/SUU;      % simple approximation
    K1 = (beta*Qi*(SXU10-alpha*SXU00-gamma*SU0')+Ri*(SYU-C*SXU));
    K2 = (Qi*beta.^2*SUU00+Ri*SUU);
    dLdKe = 2*(K1 - Ke*K2);

    %% dL/dQ
    %dLdQ = - L11 / Q^2  +  (N-1)/Q;
    dLdQ = Qi*(L1  -  (N-1));
    
    %% dL/dR 
    % dLdR = -L33 / R^2 + N/R;
    dLdR = Ri*(L3  -  N);
 

% %     %% update initial condition
% %     dLdXo = 2*(-Xs(:,1) + Xo)/Po;
% %     dLdPo = (-L5 + 1)/Po;

    dLL = 0.5*[dLdABC dLdKe dLdQ dLdR]';
% %     dLL = 0.5*[dLdABC dLdKe dLdQ dLdR dLdXo dLdPo];
end


return;





%% old implementation
C=1;

S1 = sum(Xs(2:end));
S0 = S1+Xs(1)-Xs(end);
SU0 = sum(U(:,1:end-1)')';

P = sum(Ps,3);
P11 = P - Ps(:,:,1);
P10 = sum(Pcs(:,:,2:N),3);
P00 = P - Ps(:,:,end);

SP = Xs*Xs' + P;
SP11 = SP - (Xs(:,1)*Xs(:,1)' + Ps(:,:,1));
SP10 = Xs(:,2:end)*Xs(:,1:end-1)' + sum(Pcs(:,:,2:N),3);
SP00 = Xs(:,1:end-1)*Xs(:,1:end-1)' + sum(Ps(:,:,1:N-1),3);

SXU = Xs*U';
SXU00 = SXU - Xs(:,end)*U(:,end)';
SXU10 = Xs(:,2:end)*U(:,1:end-1)';
SUU = U*U';
SUU00 = SUU - U(:,end)*U(:,end)';

SYU = Y*U';
SYY = Y*Y';
SYX = Y*Xs';



%Q and R are scalar ...
%E[(Xs(2:end)-alpha*Xs(1:end-1)-beta*Ke*U(:,1:end-1)-gamma)'Qinv(Xs(2:end)-alpha*Xs(1:end-1)-beta*Ke*U(:,1:end-1)-gamma)];
%ABC = [SP10 SXU10*Ke' S1]/[SP00 SXU00*Ke' S0; Ke*SXU00' Ke*SUU00*Ke' Ke*SU0; S0 SU0'*Ke' N-1];
ABC = [alpha beta gamma];
t1 = (SP11-ABC*[SP10'; Ke*SXU10'; S1])/Q;
t2 = (N-1)*log(abs(Q));
%E[(Y-Xs-Ke*U)'Rinv(Y-Xs-Ke*U)];
t3= (SYY + C*SP*C' + Ke*SUU*Ke'- Ke*SYU' - SYU*Ke' - C*SYX' -SYX*C' + C*SXU*Ke' + Ke*SXU'*C')/R;
t4 = N*log(abs(R));
t5 = (Xs(1)-Xo).^2/Po;
t6 = log(abs(Po));
%LL = -.5 * (t1*t1'/Q + (N-1)*log(abs(Q)) + t2*t2'/R + N*log(abs(R)) + (Xs(1)-Xo).^2/Po + log(abs(Po)));
%LL = .5 * (t1*t1'/Q + (N-1)*log(abs(Q)) + t2*t2'/R + N*log(abs(R)) + (Xs(1)-Xo).^2/Po + log(abs(Po)));
LL = .5*(t1+t2+t3+t4+t5+t6);

if isnan(LL)
    keyboard
end

if nargout > 1      % compute gradient
    %% now let's calc gradient

    %% calculate derivatives 
    %Ke = (SYU-C*SXU)/SUU;      % simple approximation
    Qi = 1/Q;   % inv(Q);
    Ri = 1/Q;   % inv(R);
    dLdABC = Qi*([SP10 SXU10*Ke' S1]-[alpha beta gamma]*[SP00 SXU00*Ke' S0; Ke*SXU00' Ke*SUU00*Ke' Ke*SU0; S0 SU0'*Ke' N-1]);
    dLdKe = (beta*Qi*(SXU10-alpha*SXU00-gamma*SU0')+Ri*(SYU-C*SXU)) - Ke*(Qi*beta.^2*SUU00+Ri*SUU);

    dLdQ = 0.5*((SP11-ABC*[SP10'; Ke*SXU10'; S1])*Qi-(N-1))*Qi;

    dLdR = 0.5*((SYY + C*SP*C' + Ke*SUU*Ke'- Ke*SYU' - SYU*Ke' - C*SYX' -SYX*C'  ...
                 + C*SXU*Ke' + Ke*SXU'*C')*Ri -N)*Ri;


    %% update initial condition
    dLdXo = (Xs(:,1)-Xo)/Ps(:,:,1);
    dLdPo = 0.5*((Xs(:,1)-Xo)*(Xs(:,1)-Xo)'/Po-1)/Po;

    dLL = [dLdABC dLdKe dLdQ dLdR dLdXo dLdPo];
    dLL = -dLL;
end
