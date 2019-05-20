function [Xp Pp Xf Pf Kf LL] = kalman_filt_1d(Y, U, A, B, C, D, Q, R, Xo, Po)   %#eml
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kalman filtering for the folling model
% X(k+1) = A*X(k) + B*U(k) + v(k)
% Y(k) = C*X(k) + D*U(k) + w(k)
% Ev = 0, Evv' = Q
% Ew = 0, Eww' = R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% to generate mex/c 
% emlmex kalman_filt_1d.m -eg {emlcoder.egs(0,[1, T]) emlcoder.egs(0,[1,T]) 0 0 0 0 0 0 0 0}
% emlmex kalman_filt_1d.m -eg {emlcoder.egs(0,[1, Inf]) emlcoder.egs(0,[Inf,1]) 0 0 0 0 0 0 0 0} -s comp_cfg


% make sure N_st == 1


N = size (Y,2);
%[N_st] = size(A,1);
% make sure N_st == 1
Xp = zeros (1,N);
Pp = zeros (1,N);
Xf = zeros (1,N);
Pf = zeros (1,N);
Kf = 0;

if (isempty(B))
    B=0;
end

LL = 0;


for i=1:N
    if i==1 % initialize
        Xp(1)=Xo;
        Pp(1)=Po;
    else
        Xp(i) = A*Xf(i-1) + B*U(:,i-1);
        Pp(i) = A*Pf(i-1)*A'+Q;
    end
    Rei=C*Pp(i)*C'+R;
    %Rei=(Rei'+Rei)/2;
    Kf=Pp(i)*C'/Rei; 
    innov=Y(i)-C*Xp(i)-D*U(:,i);
    Xf(i)=Xp(i)+Kf*innov;
    Pf(i)=Pp(i)-Kf*C*Pp(i);
    %Pf(:,:,i)=max(Pp(:,:,i)-Kf*C*Pp(:,:,i),eps);

    if (nargout >5)
        %LL= LL + log(det(Rei)) + innov'*ReiInv*innov;  % for speed up
        LL= LL + log(abs(det(Rei))) + innov'/Rei*innov;  
    end
end
LL = -.5*LL;
