function [LL] = calcNegLL(Xs, Y, U, Theta)
%function [LL dLL] = calcNegLL(Xs, Y, U, Theta)

N = length(Y);                         
M=size(U,1);
[alpha beta gamma Ke Q R Xo Po] = deal(Theta(1), Theta(2), Theta(3), ...
                                       Theta(4:M+3), Theta(M+4), Theta(M+5), ...
                                       Theta(M+6), Theta(M+7));   
Ke = Ke';
t1 = Xs(2:end) - alpha*Xs(1:end-1) - beta*Ke*U(:,1:end-1) - gamma;
t2 = Y-Xs-Ke*U;

%Q and R are scalar ...
LL = .5 * (t1*t1'/Q + (N-1)*log(abs(Q)) + t2*t2'/R + N*log(abs(R)) + Xo.^2/Po + log(abs(Po)));
%+ Xs(1).^2/;

if isnan(LL)
    keyboard
end

%% let's calc derivatives
% 
% return;
% 
% 
% %% Calc needed values 
% % for updating alpha and beta
% S1 = sum(Xs(2:end));
% S0 = S1+Xs(1)-Xs(end);
% SU0 = sum(U(:,1:end-1)')';
% 
% S11 = 0;    %Xs(:,1)*Xs(:,1)' + Ps(:,:,1);
% S10 = 0;    %Xs(:,1)*Xo' + Pcs(:,:,1);
% S00 = 0;    %Xo*Xo' + Po;
% SXU10 = 0;
% SXU00 = 0;
% SUU00 = 0;
% 
% % for updating Ke
% SYU = Y(:,1)*U(:,1)';
% 
% % for update R
% SYY = Y(:,1)*Y(:,1)';% + C*Ps(:,:,i)*C';
% SYX = Y(:,1)*Xs(:,1)';
% 
% for i=2:N
%     % for updating alpha and beta
%     S11 = S11 + Xs(:,i)*Xs(:,i)' ;%+ Ps(:,:,i);
%     S10 = S10 + Xs(:,i)*Xs(:,i-1)' ;%+ Pcs(:,:,i);
%     S00 = S00 + Xs(:,i-1)*Xs(:,i-1)' ;%+ Ps(:,:,i-1);
% 
%     SXU10 = SXU10 + Xs(:,i)*U(:,i-1)';
%     SXU00 = SXU00 + Xs(:,i-1)*U(:,i-1)';
%     SUU00 = SUU00 + U(:,i-1)*U(:,i-1)';
% 
%     % for updating Ke
%     SYU = SYU + Y(:,i)*U(:,i)';
% 
%     % for update R
%     SYY = SYY + Y(:,i)*Y(:,i)';% + C*Ps(:,:,i)*C';
%     SYX = SYX + Y(:,i)*Xs(:,i)';
% 
% end
% 
% 
% %% calc derivatives 
% % 1) dLL/d[alpha beta gamma]
% ABC=[alpha beta gamma];
% Qi = 1/Q;   %inv(Q);
% Ri = 1/R;   %inv(R);
% dLLdABC = -Qi*([S10 SXU10*Ke' S1] - ABC*[S00 SXU00*Ke' S0; Ke*SXU00' Ke*SUU00*Ke' Ke*SU0; S0 SU0'*Ke' N-1]);
% 
% % 2)dLL/dKe
% SXU = SXU00 + Xs(:,end)*U(:,end)';
% SUU = (SUU00 + U(:,end)*U(:,end)');
% %Ke = (SYU-C*SXU)/SUU;      % simple approximation
% dLLdKe = -(beta*Qi*(SXU10-alpha*SXU00-gamma*SU0')+Ri*(SYU-SXU)) ...
%         + Ke*(Qi*beta.^2*SUU00+Ri*SUU);
% 
% % 3) dLL/dQ
% dLLdQ=(S11-ABC*[S10'; Ke*SXU10'; S1])*Qi.^2-(N-1)*Qi;  
% 
% % 4) dLL/dR
% S = (S11+Xs(:,1)*Xs(:,1)')   ;%+ Ps(:,:,1));
% SXU = SXU00 + Xs(:,end)*U(:,end)';
% SUU = SUU00 + U(:,end)*U(:,end)';
% dLLdR = (SYY + S + Ke*SUU*Ke'- Ke*SYU' - SYU*Ke' - SYX' -SYX  ...
%      + SXU*Ke' + Ke*SXU')*Ri.^2 - N*Ri;
% 
% % 5) dLL/dXo
% dLLdXo = -(Xs(:,1)-Xo)/Po;
% 
% % 6) dLL/dPo
% dLLdPo = dLLdXo.^2 -1/Po;
% 
% 
% dLL = [dLLdABC dLLdKe dLLdQ dLLdR dLLdXo dLLdPo];
