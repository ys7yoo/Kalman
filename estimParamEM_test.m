clear

load Y.txt
load I.txt

paramLow =[0.9, 1e-6, -0.1, 0.1, 0.1];
paramHigh =[0.99, 0.01, 0.1, 0.5, 0.5];
initPrior = 0.01;
    
nke =6
[paramEst KeEst LL LLs initParam initKe allParam allKe] = estimParamEM(Y, I, paramLow, paramHigh, initPrior, nke, 200, 8)
LLtrials = max(LLs);                % ending LL at each trial
itrTrials = sum(~isnan(LLs));       % number of iterations per each trial
%%
numRow=3;
clf
subplot(numRow,1,1)
plot(LLs)
set(gca,'xlim', [1 max(itrTrials)])
box off

subplot(numRow,2,3)
%hist(LLtrials)
plot (LLtrials, 'o--')
title ('final LLs')
xlabel ('trial')
box off

subplot(numRow,2,4)
plot (itrTrials, 'o--')
title ('iterations')
%hist(itrTrials)
%title ('histogram of iteration before convergence')
xlabel ('trial')
box off
% % Results are like the following
% paramEst =
%     0.9911
%     0.0082
%    -0.4243
%     0.1075
%     0.0975
% 
% KeEst =
%    -0.0003    0.7998    0.4016    0.2008    0.1017    0.0512    

subplot(numRow,2,5)
if exist('initAt','var')
    plot(initAt','g+:'); hold on
    plot(allAt','bo:');

    title ('initial and estimated aT')
else
    plot(initKe','g+:'); hold on
    plot(allKe','bo:');

    title ('initial and estimated Ke')
end

xlabel('trial')
box off

subplot(numRow,2,6)
plot(initParam','g+:'); hold on
plot(allParam','bo:'); hold on
xlabel('trial')
title ('initial and estimated param')
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADD the same routine for multi-tap fitting 
clear
load multiY.txt
load multiU.txt     % input is already filtered!
Y = multiY;
U = multiU;
% Y = Y';
% I = I';
nke = 6;

dimX = 2;

AtLow = 0.9*ones(1,dimX);
AtHigh = 1.1*ones(1,dimX);

ParamLow = [1e-6, -0.1, 0.1, 0.1];
ParamHigh = [0.01, 0.1, 0.5, 0.5];

initPrior = 10;

[AtEst, paramEstMulti, LL, LLs, initAt, initParam, allAt, allParam] = estimParamEM_multi(Y, U, AtLow, AtHigh, ParamLow, ParamHigh, initPrior, 20, 5)
LLtrials = max(LLs);                % ending LL at each trial
itrTrials = sum(~isnan(LLs));       % number of iterations per each trial


%% Plot initial and estimated parameters 
