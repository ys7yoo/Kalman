function [Xs Ps Pcs Xf Pf Xp Pp LL] = kalman_smth(Y, U, A, B, C, D, Q, R, Xo, Po) %#eml
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kalman smoothing for the following model
% X(k+1) = A*X(k) + B*U(k) + v(k)
% Y(k) = C*X(k) + D*U(k) + w(k)
% Ev = 0, Evv' = Q
% Ew = 0, Eww' = R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(Y,2);
[N_st] = size(A,1);
Xs = zeros (N_st,N);
Ps = zeros (N_st,N_st, N);
J = zeros (N_st,N_st, N-1);

% call filter first
if nargin > 7
    [Xp Pp Xf Pf Kn LL] = kalman_filt(Y, U, A, B, C, D, Q, R, Xo, Po);
else
    [Xp Pp Xf Pf Kn] = kalman_filt(Y, U, A, B, C, D, Q, R, Xo, Po);
end
    

% initialize at the end
Xs(:,N) = Xf(:,N);
Ps(:,:,N) = Pf(:,:,N);
for k=N:-1:2
    %J(:,:,k-1)=Pf(:,:,k-1)*A'*inv(Pp(:,:,k)); % for speed up
    J(:,:,k-1)=Pf(:,:,k-1)*A'/Pp(:,:,k); % for speed up
    Xs(:,k-1)= Xf(:,k-1) + J(:,:,k-1)*(Xs(:,k)-Xp(:,k));
    Ps(:,:,k-1)= Pf(:,:,k-1) + J(:,:,k-1)*(Ps(:,:,k)-Pp(:,:,k))*J(:,:,k-1)';
%     if (Ps(1,1,k-1) <0)
%         keyboard
%     end
end

% Lag-One Covariance Smoothers (Pcs=P_{t,t-1}^n)
Pcs = zeros (N_st,N_st, N);
Pcs(:,:,N) = (eye(N_st)-Kn*C)*A*Pf(:,:,N-1);

for k=N:-1:3
   Pcs(:,:,k-1)=Pf(:,:,k-1)*J(:,:,k-2)' + ...
                J(:,:,k-1)*(Pcs(:,:,k)-A*Pf(:,:,k-1))*J(:,:,k-2)';
end
Pcs(:,:,1) = NaN;
           