function [LL, dLL] = calcNegLLparamLog(Xs, Y, U, Ps, Pcs, Xo, Po, Theta)

% The goal is to calculate negative log likelihood, -log P(Y|Theta)
% given Y, Theta, and P(X|Y,Theta) (Xs, Ps, Pcs), 
% In this implementation, Q and R is parameterized as Q = exp(-q) and R = exp(-r)


N = length(Y);                         
M=size(U,1);

[alpha beta gamma Ke q r] = deal(Theta(1), Theta(2), Theta(3), ...
                                       Theta(4:M+3), Theta(M+4), Theta(M+5));
Q = exp(-q);
R = exp(-r);
                                   
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
L2 = (N-1)*q;
%L2 = (N-1)*log(abs(Q));

L33 = (Y-Xs-KeU)*(Y-Xs-KeU)' + P;
L3 = L33/R;
%L4 = N*log(abs(R));
L4 = N*r;

L5 = (Xs(1)-Xo).^2/Po;
L6 = log(abs(Po));

LL = .5*(L1+L2+L3+L4+L5+L6);


if nargout >1 % calc gradient too
    Qi = exp(q);
    Ri = exp(r);
    %Qi = 1/Q;
    %Ri = 1/R;
    
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
    dLdABC = 2*Qi*( -A1 + [alpha beta gamma]*A2);
    
    
    %% calculate new Ke
    %Ke = (SYU-C*SXU)/SUU;      % simple approximation
    K1 = (beta*Qi*(SXU10-alpha*SXU00-gamma*SU0')+Ri*(SYU-C*SXU));
    K2 = (Qi*beta.^2*SUU00+Ri*SUU);
    dLdKe = 2*(-K1 + Ke*K2);

    %% dL/dQ
    %dLdQ = - L11 / Q^2  +  (N-1)/Q;
    %dLdQ = Qi*(- L1  +  (N-1));
    dLdq = Qi*(- L11)  +  (N-1);
    
    %% dL/dR 
    % dLdR = -L33 / R^2 + N/R;
    %dLdR = Ri*(- L3  +  N);
    dLdr = Ri*(- L3)  +  N;


    dLL = 0.5*[dLdABC dLdKe dLdq dLdr]';
end


return;




