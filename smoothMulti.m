function [Vest] = smoothMulti(At, paramEstMulti, initPrior, Y, U)


% paramEstMulti = [beta, Vo, q, R,]
%smoothMulti(AtEst, paramEstMulti(1), 
%% Smoothe with estimated parameters and calc MSE
% the last input V is the  true membrane potential if accessible
% 
% AtEst for At
% paramEstMulti(1) for beta 
% paramEstMulti(2) for Vo
% paramEstMulti(3); for q
% paramEstMulti(4); for R
%% Smoothe with estimated parameters and calc MSE

dimX = length(At);

% smooth here     
A = zeros(dimX,dimX);
A(1,:) = At;
A(2:end,1:end-1) = eye(dimX-1);

%                    B = zeros(dimX,1);
%                     B(1,1) = paramEstMulti(1);      % beta
B = paramEstMulti(1);

Vo = paramEstMulti(2);

C = [1 zeros(1,dimX-1)];
D = 1;

Q = zeros(dimX,dimX);
Q(1,1) = paramEstMulti(3);

R = paramEstMulti(4);

Xo = zeros(dimX,1);

if numel(initPrior)==1
    initPrior = initPrior*eye(dimX);
end

[Xs] = kalman_smth(Y'-Vo, U', A, B, C, D, Q, R, Xo, initPrior);
%[Xs Ps Pcs] = kalman_smth(Y'-Vo, U', A, B, C, D, Q, R, Xo, initPrior);
Vest = Xs(1,:)' + Vo;

