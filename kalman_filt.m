function [Xp Pp Xf Pf Kf LL] = kalman_filt(Y, U, A, B, C, D, Q, R, Xo, Po) %#eml
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kalman filtering for the folling model
% X(k+1) = A*X(k) + B*U(k) + v(k)
% Y(k) = C*X(k) + D*U(k) + w(k)
% Ev = 0, Evv' = Q
% Ew = 0, Eww' = R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size (Y,2);
[N_ms N_st] = size(C);
Xp = zeros (N_st,N);
Pp = zeros (N_st,N_st, N);
Xf = zeros (N_st,N);
Pf = zeros (N_st,N_st, N);
Kf = zeros (N_st,N_ms);
if (isempty(B))
    B=0;
end

LL = 0;

for i=1:N
    if i==1 % initialize
        Xp(:,1)=Xo;
        Pp(:,:,1)=Po;
    else
        Xp(:,i) = A*Xf(:,i-1) + B*U(:,i-1);
        Pp(:,:,i) = A*Pf(:,:,i-1)*A'+Q;
    end
    Rei=C*Pp(:,:,i)*C'+R;
    Rei=(Rei'+Rei)/2;
    %ReiInv=inv(Rei);
    %Kf=Pp(:,:,i)*C'*ReiInv;    % for speed up
    Kf=Pp(:,:,i)*C'/Rei; 
    innov=Y(:,i)-C*Xp(:,i)-D*U(:,i);
    Xf(:,i)=Xp(:,i)+Kf*innov;
    Pf(:,:,i)=Pp(:,:,i)-Kf*C*Pp(:,:,i);
    %Pf(:,:,i)=max(Pp(:,:,i)-Kf*C*Pp(:,:,i),eps);

%     if (Pf(1,1,i) <0)
%         keyboard
%     end

    if (nargout >5)
        %LL= LL + log(det(Rei)) + innov'*ReiInv*innov;  % for speed up
        LL= LL + log(abs(det(Rei))) + innov'/Rei*innov;  
    end
end

LL = -.5*LL;
