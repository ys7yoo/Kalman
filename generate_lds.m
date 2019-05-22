function [X,Y] = generate_lds(U, A, B, C, D, Q, R, initX, initV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate sample trajectory of LDS:
% X(k+1) = A*X(k) + B*U(k) + v(k)
% Y(k) = C*X(k) + D*U(k) + w(k)
% Ev = 0, Evv' = Q
% Ew = 0, Eww' = R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = size(U,2);
[os, ss] = size(C);
X = zeros(ss, T);
Y = zeros(os, T);

sQ = sqrt(diag(Q));
sR = sqrt(diag(R));

% draw initial state from initX and initV and measure it
X(:,1) = initX + sqrt(initV)*randn(ss,1);  % iniV must be diagonal

if isempty(U)
    Y(:,1) = C*X(:,1) + randn(os,1).*sR;
else
    Y(:,1) = C*X(:,1) + D*U(:,1) + randn(os,1).*sR;
end

% Another way of initialization
% draw initial state
%Xo = initX + randn(ss,1).*sqrt(diag(initV));
% X(:,1) = A*Xo + random('Normal', 0, diag(Q));
% Y(:,1) = C*X(:,1)  + random('Normal', 0, diag(R));
if isempty(U)
    for t=2:T
      X(:,t) = A*X(:,t-1) + sQ.*randn(ss,1);
      Y(:,t) = C*X(:,t)  + sR.*randn(os,1);
    end
else
    for t=2:T
      X(:,t) = A*X(:,t-1) + B*U(:,t-1) + sQ.*randn(ss,1);
      Y(:,t) = C*X(:,t)  + D*U(:,t) + sR.*randn(os,1);
    end
end

