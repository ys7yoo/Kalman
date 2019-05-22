%% Choose param of a linear dynamic system 

% dynamics
A = 0.9*randn(2,2)
B = zeros(2,1)
Q = 0.01*eye(2,2)

% measurements
C = randn(1,2)
D = randn(1,1)
R = 0.01

% control input
U = zeros(1,40)

%% generate (X,Y) according to the linear dynamics system
initX = zeros(2,1)
initV = 0.0001*eye(2,2)
[X,Y] = generate_lds(U, A, B, C, D, Q, R, initX, initV)

figure(1)
clf
subplot(211)
plot(X')
xlabel('n')
ylabel('X_n')

subplot(212)
plot(Y')
xlabel('n')
ylabel('Y_n')



%% estimate X from Y using Kalman filter
Xo = initX;
Po = initV;
[Xp Pp Xf Pf Kf LL] = kalman_filt(Y, U, A, B, C, D, Q, R, Xo, Po)


figure(2)
clf
subplot(211)
plot(Xf(1,:)'); hold on
plot(X(1,:)', '--')
ylabel('X_n(1)')

legend('estim', 'true')

subplot(212)
plot(Xf(2,:)'); hold on
plot(X(2,:)', '--')
ylabel('X_n(2)')

figure(3)
clf
plot(X(:), Xf(:), 'o')
xlabel('true X_n')
ylabel('estimated X_n')

axis equal


%% estimate param from (X,Y)
param_init.A = A+0.1*randn(2,2);
param_init.B = B+0.1*randn(2,1);
param_init.C = C+0.1*randn(1,2);
param_init.D = D;
param_init.Q = Q;
param_init.R = R;
param_init.Xo = Xo;
param_init.Po = Po;


[param_estim] = em_kalman_abcd(Y, U, param_init)

