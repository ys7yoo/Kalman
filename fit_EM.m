function [paramEst paramInit Klin ] = fit_EM(Iapp, Y, initGuess, FIT, paramTrue)

if nargin <5
    paramTrue = [];
end

%% T = length(Iapp);

%% Plot input and outputs
% tt = (1:T)*dt;
% idx = 1:(min(1000,T));
% clf;
% subplot(311);
% plot(tt(idx),Iapp(idx),'r');hold on
% plot(tt(idx),Iinj(idx));
% legend('I_{app}','I_{inj}');
% xlabel('ms')
% ylabel('uA/cm^2')
% 
% 
% subplot(312);
% plot(tt(idx),V(idx));
% title('V_m');
% xlabel('ms')
% ylabel('mV')
% 
% subplot(313);
% plot(tt(idx),[m(idx) n(idx) h(idx)]);
% %plot(tt,[m n h]);
% title('Gating variables');
% xlabel('ms');
% ylabel('1/ms')
% legend ('m', 'n', 'h')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate alpha,beta,Ke,Q,R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial guess of kernel using simple linear regression
%% changed to consider mean (May 24, 2011)
II=[stackCols(Iapp,FIT.M,0); ones(1,size(Iapp,1)) ] ;
KlinVo = II'\Y(:);
Klin = KlinVo(1:end-1)';
Vo = KlinVo(end);

%II=stackCols(Iapp,FIT.M,0);
%Klin = II'\Y(:);
%Klin = Klin' - Klin(1); % simple approximation  % not good for real data
%Klin = Klin' - max(0,min(Klin));  %remove large positive bias
%Klin = Klin' - min(Klin);  %remove large positive bias
%Klin = Klin';

% % Klin = makeStimRows(Iapp,FIT.M)\Y';
% % Klin = flipud(Klin)';
% % Klin = Klin - Klin(1);
    

% plot ([Klin' paramTrue.Ke'])


clf; 
clear paramInit paramEstim LLtrace;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for trial = 1:FIT.num_trial
    % initialize
    if (isempty(initGuess))
        paramInit(trial).alpha = .9 + .1*rand();
        paramInit(trial).beta = .1/5.55 + 0.01*rand();
        paramInit(trial).gamma = 0.3*0.1*10.6 + 0.3*rand();
        paramInit(trial).C = 1;
        paramInit(trial).Ke = Klin + 0.1*randn(1,FIT.M);
        %paramInit(trial).Ke = [0 initGuess.sumKe/(FIT.M-2) + randn(1,FIT.M-2) 0];
    %     paramInit(trial).Q = Qinit*(.5+rand()); %min(0.25*rand+.25,Rinit);
    %     paramInit(trial).R = Rinit*(.5+rand());
        paramInit(trial).Q = initGuess.Q;
        paramInit(trial).R = initGuess.R
        paramInit(trial).Xo = Vo; 
        paramInit(trial).Po = 5;
        
        % %     % start some of the value from the original value (for debug)
% %     paramInit(trial).alpha = 0.97;
% %     paramInit(trial).beta = .1/5.55;
% %     paramInit(trial).gamma = 0.3*0.1*10.6;
% %     paramInit(trial).C = 1;    
% %     paramInit(trial).Ke =[0 0 0.4 3.5 1.5 0.15 0];
% %     paramInit(trial).Q = 1;
% %     paramInit(trial).R = 1;
% %     paramInit(trial).Xo = -60;
% %     paramInit(trial).Po = 1;
    else
% %         switch initGuess.type
% %             case 'fixed'
% %                 paramInit(trial) = initGuess;
% %                 paramInit(trial).C = 1;
% %                 
% %             otherwise   %  add perturbation to initial guess
% %                 
                
                % change initial valu with range (Mar. 13, 2011)
                if (numel(initGuess.alpha)==1)
                    paramInit(trial).alpha = random('unif', initGuess.alpha(1),1);
                else
                    paramInit(trial).alpha = random('unif', initGuess.alpha(1),initGuess.alpha(2));
                end
                if exist('FIT.fixBeta','var') & FIT.fixBeta
                    paramInit(trial).beta = initGuess.beta;
                else
                    if (numel(initGuess.beta)==1)
                        paramInit(trial).beta = initGuess.beta + 0.01*randn();
                    else
                        paramInit(trial).beta = random('unif', initGuess.beta(1),initGuess.beta(2));
                    end
                end
                if (numel(initGuess.gamma)==1)      %Vo + sqrt(initGuess.Po)*randn();   %% not a goo idea!   
                    paramInit(trial).gamma = initGuess.gamma + 0.1*randn();
                else
                    paramInit(trial).gamma = random('unif', initGuess.gamma(1),initGuess.gamma(2));
                    
                end

        %         paramInit(trial).alpha = initGuess.alpha + (1-initGuess.alpha)*rand();
        %         paramInit(trial).beta = initGuess.beta*(5*rand());
        %         paramInit(trial).gamma = initGuess.gamma*(3*rand());
                paramInit(trial).C = 1;
                paramInit(trial).Ke = max(Klin + 0.1*randn(1,FIT.M),0);
                if (numel(initGuess.Q)==1)      
                    paramInit(trial).Q = initGuess.Q;
                else
                    paramInit(trial).Q = random('unif', initGuess.Q(1),initGuess.Q(2));
                end
                if (numel(initGuess.R)==1)      
                    paramInit(trial).R = initGuess.R;
                else
                    paramInit(trial).R = random('unif', initGuess.R(1),initGuess.R(2));
                end
% %                 paramInit(trial).Q = initGuess.Q;
% %                 paramInit(trial).R = initGuess.R;

                
                % prior
                paramInit(trial).Xo = initGuess.Xo;
                paramInit(trial).Po = initGuess.Po;

    % %         % initGuess without noise (for debug)
    % %         paramInit(trial).alpha = initGuess.alpha;
    % %         paramInit(trial).beta = initGuess.beta;
    % %         paramInit(trial).gamma = initGuess.gamma;
    % %         paramInit(trial).Q = initGuess.Q;
    % %         paramInit(trial).R = initGuess.R;
    % %         paramInit(trial).Xo = initGuess.Xo;
    % %         paramInit(trial).Po = initGuess.Po;
    % % 
    % %         paramInit(trial).C = 1;
    % %         paramInit(trial).Ke = Klin;
% %         end

           
            

    end    
end

LLend = zeros(FIT.num_trial,1);
%XXs = zeros(FIT.num_trial,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run FIT FIT.num_trial tiems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ticID = tic;
%parfor trial = 1:FIT.num_trial
for trial = 1:FIT.num_trial
    % estimate parameters
    tic;
    [paramEstim(trial) LL] = em_kalman(Y, Iapp, paramInit(trial), FIT, trial, paramTrue);        % use Iapp instead of II because I can easily generate Ke*Iapp with the filter command
    
    LLtrace(trial).LL = LL;   
    LLend(trial) = LL(end);
    fprintf('%dth tiral: LL=%.2e (%d itr, %.1fs)\n', trial, LL(end), length(LL),toc)
    
    

end
fprintf('\nTotal %d trials: LLmax=%.2e (%.1fs)\n', FIT.num_trial, max(LLend), toc(ticID))

%% plot sammary accross trials

%% draw
clf;
for trial = 1:FIT.num_trial
    subplot(331); hold on;
    semilogy(real(LLtrace(trial).LL));
    %semilogy(real([LL negLL]));
    %title ('Log likelihood')
    title ('L(\Theta)')
    xlabel ('EM iteration')
    
% %     subplot(253); hold on;
% %     semilogy(diff(real(LLtrace(trial).LL)));
% %     %semilogy(real([LL negLL]));
% %     %title ('Log likelihood')
% %     title ('\Delta L(\Theta)')
% %     xlabel ('EM iteration')
    
    subplot(332); hold on;
    plot(trial, LLtrace(trial).LL(1), 'x')
    plot(trial, LLend(trial),'o');
    %legend ('initia', 'final')
    title ('Final Log likelihood (x: initial, o: final)')
    xlabel ('Trial')
    
    
    
% %     subplot(254); hold on; 
% %     plot(trial, paramEstim(trial).gamma, 'o'); title ('\gamma'); xlabel('trial')
% % 
% %     
end


% plot initial param
%paramInit(trial).alpha

subplot(334);
plot([paramInit(:).alpha], 'xk'); hold on
subplot(335);
plot([paramInit(:).beta], 'xk'); hold on
subplot(336);
plot([paramInit(:).gamma], 'xk'); hold on

subplot(338);
plot([paramInit(:).Q], 'xk'); hold on
subplot(339);
plot([paramInit(:).R], 'xk'); hold on

% estimated params
subplot(334);
plot([paramEstim(:).alpha], 'o:'); title ('\alpha'); xlabel('trial');box off
subplot(335);
plot([paramEstim(:).beta], 'o:'); title ('\beta'); xlabel('trial'); box off
subplot(336);
plot([paramEstim(:).gamma], 'o:'); title ('\gamma'); xlabel('trial'); box off
subplot(337);
KeHat = reshape([paramEstim(:).Ke], FIT.M,[]);
%plot(reshape([paramEstim(:).Ke], FIT.M,[]),'o-')
if (size(KeHat,2)==1)
    plot(KeHat, 'o:'); box off
else
    errorbar(mean(KeHat'), std(KeHat')); box off
end
xlabel('trial'); box off
subplot(338);
plot([paramEstim(:).Q], 'o:'); title ('Q'); xlabel('trial'); box off
subplot(339);
plot([paramEstim(:).R], 'o:'); title ('R'); xlabel('trial'); box off
% subplot(259);
% plot([paramEstim(:).Xo], 'o:'); title ('Xo'); xlabel('trial'); box off
% subplot(2,5,10);
% plot([paramEstim(:).Po], 'o:'); title ('Po'); xlabel('trial'); box off




% plot true param
if (~isempty(paramTrue))
    subplot(334);  hold on
    plot(repmat(paramTrue.alpha,FIT.num_trial), '+r--'); 
    subplot(335); hold on
    plot(repmat(paramTrue.beta,FIT.num_trial), '+r--'); 
    subplot(336); hold on
    plot(repmat(paramTrue.gamma,FIT.num_trial), '+r--'); 
    subplot(337); hold on
    plot(paramTrue.Ke, '+r--'); 
    subplot(338);  hold on
    plot(repmat(paramTrue.Q,FIT.num_trial), '+r--'); 
    set(gca,'ylim', [0 max(paramTrue.Q*2, max([paramEstim(:).Q]))]);
    subplot(339); hold on
    plot(repmat(paramTrue.R,FIT.num_trial), '+r--'); 
    set(gca,'ylim', [0 max(paramTrue.R*2, max([paramEstim(:).R]))]);
% %     subplot(259); hold on
% %     plot(repmat(paramTrue.Xo,FIT.num_trial), '+r--'); 
% %     subplot(2,5,10); hold on
% %     plot(repmat(paramTrue.Po,FIT.num_trial), '+r--'); 
% %     set(gca,'ylim', [0 max(paramTrue.Po*2, max([paramEstim(:).Po]))]);
end

set(gcf, 'paperposition', [0 0 15 6]) 
set(gcf, 'papersize', [15 6]) 
saveas(1, sprintf('params_EM.pdf'));

%% choose parameter set with maximum LL
[~, idx] = max(LLend);
paramEst = paramEstim(idx);
%paramEst.Xs = XXs(idx,:);


%% now filter/smoother with the best param
[paramEst.Xs paramEst.Ps , ~, ~, paramEst.Xf paramEst.Pf] = kalman_smth_1d(Y-paramEst.gamma, filter(paramEst.Ke,1,Iapp)', paramEst.alpha, paramEst.beta, 1, 1, paramEst.Q, paramEst.R, initGuess.Xo, initGuess.Po);



return



% below: analyze result 

Ifil = filter(Ke, 1, Iinj);
IfilEst = filter(paramEst.Ke, 1, Iinj);


% original parameters
gL = .1;
Cm =1;
Vo = -65;

alpha = 1-dt*gL/Cm;
beta = gL/Cm/sum(Ke);
gamma = gL*dt*Vo/Cm;



% save result
filename = sprintf('estimation_nn_R%.0f_itr%d.mat', R, FIT.max_iter);
save(filename)



%% draw figures
subplot(612);hold off
plot (dt*(1:T), Iinj,'k--'); hold on;
plot (dt*(1:T), Ifil,'r');
plot (dt*(1:T), IfilEst);
%plot (dt*(1:T),[Iinj(:) Ifil(:)])
legend ('I_{app}', 'True I_{inj}', 'Estimated I_{inj}')
xlabel('ms')
ylabel('\muA')
title ('Input current')
subplot(613);hold off
plot (dt*(1:T),paramEst.Xs(1,:),'r'); hold on
plot (dt*(1:T),X(1,:))
legend ('True V_m', 'Estimated V_m')
xlabel('ms')
ylabel('mV')
title (sprintf('Membrane voltage'))
subplot(614)
plot (dt*(1:T),Ke*II,'r'); hold on
plot (dt*(1:T),paramEst.Ke*II)
legend ('True V_e', 'Estimated V_e')
xlabel('ms')
ylabel('mV')
title (sprintf('Electrode voltage'))
%
% LDS parameters
subplot(6,4,17);hold off
plot (alpha, '+r'); hold on
plot (paramEst.alpha, 'o')
%legend ('True', 'Estimated')
title ('\alpha')
set(gca,'ylim', [min([0; paramEst.alpha]) 2*max([alpha; paramEst.alpha])]);
% compare beta
subplot(6,4,18);hold off
plot ([beta], '+r'); hold on
plot ([paramEst.beta], 'o')
%legend ('True', 'Estimated')
title ('\beta')
set(gca,'ylim', [min([0; paramEst.beta]), 2*max([beta; paramEst.beta])]);
% compare gamma
subplot(6,4,19);hold off
plot ([gamma], '+r'); hold on
plot ([paramEst.gamma], 'o')
%legend ('True', 'Estimated')
title ('\gamma')
set(gca,'ylim', [2*min([gamma; paramEst.gamma]), max([0; paramEst.gamma]), ]);
% compare kernel
subplot(6,4,20); hold off
plot (Ke', '+-r'); hold on
%plot (Klin, 'x-');
plot (paramEst.Ke, 'o-');
%legend ('True', 'Estimated')
title ('Ke')
%
% original parameter
subplot(6,4,21);hold off
plot (Cm, '+r'); hold on
plot (paramEst.Cm, 'o')
%legend ('True', 'Estimated')
title ('C')
set(gca,'ylim', [min([0; paramEst.C]), 2*max([C; paramEst.C])]);
% compare gL
subplot(6,4,22);hold off
plot (gL, '+r'); hold on
plot (paramEst.gL, 'o')
%legend ('True', 'Estimated')
title ('g_L')
set(gca,'ylim', [min([0; paramEst.gL]) 2*max([gL; paramEst.gL])]);
% compare gamma
subplot(6,4,23);hold off
plot (Vo, '+r'); hold on
plot (paramEst.Vo, 'o')
%legend ('True', 'Estimated')
title ('Vo')
set(gca,'ylim', [2*min([Vo; paramEst.Vo]), max([0; paramEst.Vo]), ]);
% compare kernel
subplot(6,4,24);hold off
plot (Ke', '+-r'); hold on
plot (paramEst.Ke, 'o-');
legend ('True', 'Estimated')
title ('Ke')


%%
filename = sprintf('estimation_nn_R%.0f_itr%d.png',R,FIT.max_iter);
saveas(1, filename)

%%
toc

return

%% draw separate figures

%% 0) input and output
clf;
subplot(211)
plot (dt*(1:T), Iinj);
%plot (dt*(1:T),[Iinj(:) Ifil(:)])
%legend ('I_{inj}', 'I_{filterd}')
xlabel('ms')
ylabel('nA')
title ('Input current')
set(gca,'xlim',[0 100])
subplot(212)
plot (dt*(1:T),Y(:),'k'); hold on;
plot (dt*(1:T),X(1,:),'b--'); hold on
plot (dt*(1:T),(Ke*II),'r--');
legend ('Vrec', 'Vm', 'Ve')
xlabel('ms')
ylabel('mV')
title ('Measured voltage')
set(gca,'xlim',[0 100])
%title (sprintf('Membrane and electrode voltage (Q=%.0e)',Q(1,1)))
% subplot(413)
% plot (dt*(1:T),Y(:))
% xlabel('ms')
% ylabel('mV')
% title (sprintf('Measured voltage (Vm+Ve+noise)(R=%.0e)',R))

filename = sprintf('generated_Q%.0f_R%.0f_N%d_itr%d.png', Q,R,T,FIT.max_iter);
saveas(1, filename)


%% 1) draw LL traces
clf
hold on;
for trial = 1:length(LLtrace)
    semilogy(real(LLtrace(trial).LL));
end
title ('L(\Theta)')
xlabel ('Iteration')

%% 2) Parameter estimation of LDS
disp('Estimated parameters for LDS')
disp(sprintf('alpha=%.3f', paramEst.alpha))
disp(sprintf('beta=%.3f', paramEst.beta))
disp(sprintf('gamma=%.3f', paramEst.gamma))
% compare alpha
clf
subplot(141);hold off
plot (alpha, '+r'); hold on
plot (paramEst.alpha, 'o')
%legend ('True', 'Estimated')
title ('\alpha')
set(gca,'ylim', [min([0; paramEst.alpha]) 2*max([alpha; paramEst.alpha])]);
% compare beta
subplot(142);hold off
plot ([beta], '+r'); hold on
plot ([paramEst.beta], 'o')
%legend ('True', 'Estimated')
title ('\beta')
set(gca,'ylim', [min([0; paramEst.beta]), 2*max([beta; paramEst.beta])]);
% compare gamma
subplot(143);hold off
plot ([gamma], '+r'); hold on
plot ([paramEst.gamma], 'o')
%legend ('True', 'Estimated')
title ('\gamma')
set(gca,'ylim', [2*min([gamma; paramEst.gamma]), max([0; paramEst.gamma]), ]);
% compare kernel
subplot(144); hold off
plot (Ke', '+-r'); hold on
%plot (Klin, 'x-');
plot (paramEst.Ke, 'o-');
legend ('True', 'Estimated')
title ('Ke')

%% 3) Kernel only

paramEst.alpha
clf
plot (Ke', '+-r', 'MarkerSize',10); hold on
%plot (Klin, 'x-');
plot (paramEst.Ke, 'o-',  'MarkerSize',10);
legend ('True', 'Estimated')
title ('Ke')
set(gca,'ylim', [-0.1 4])

%% 4) Original parameters 
disp('Estimated parameters for LDS')
disp(sprintf('C=%.3f', paramEst.Cm))
disp(sprintf('g_L=%.3f', paramEst.gL))
disp(sprintf('V_o=%.3f', paramEst.Vo))

% original parameter
clf
subplot(141);hold off
plot (Cm, '+r'); hold on
plot (paramEst.Cm, 'o')
%legend ('True', 'Estimated')
title ('Cm')
set(gca,'ylim', [min([0; paramEst.Cm]), 2*max([Cm; paramEst.Cm])]);
% compare gL
subplot(142);hold off
plot (gL, '+r'); hold on
plot (paramEst.gL, 'o')
%legend ('True', 'Estimated')
title ('g_L')
set(gca,'ylim', [min([0; paramEst.gL]) 2*max([gL; paramEst.gL])]);
% compare gamma
subplot(143);hold off
plot (Vo, '+r'); hold on
plot (paramEst.Vo, 'o')
%legend ('True', 'Estimated')
title ('Vo')
set(gca,'ylim', [2*min([Vo; paramEst.Vo]), max([0; paramEst.Vo]), ]);
% compare kernel
subplot(144);hold off
plot (Ke', '+-r'); hold on
plot (paramEst.Ke, 'o-');
legend ('True', 'Estimated')
title ('Ke')

%% Compare estimation results


% 1) compensate by AEC
VeAEC = AEC.Keff(1:FIT.M)'*II;
VmAEC = Y-VeAEC;
mseVeAEC = (Ve-VeAEC)*(Ve-VeAEC)'/T
mseVmAEC = (X-VmAEC)*(X-VmAEC)'/T

% 2) compensage by Ke estimated by FIT
VeEMKe = paramEst.Ke*II;
VmEMKe = Y-VeEMKe;
mseVeEMKe = (Ve-VeEMKe)*(Ve-VeEMKe)'/T
mseVmEMKe = (X-VmEMKe)*(X-VmEMKe)'/T

% 3) smoothing by alpha, beta, gamma and Ke estimated by FIT
%[VmEMsm Ps] = kalman_smth(Y, IIext, paramEst.alpha,
%[paramEst.beta*paramEst.Ke paramEst.gamma], C, [paramEst.Ke 0], paramEst.Q, paramEst.R, paramEst.Xo, paramEst.Po);
mseVmEMsm = (X-paramEst.Xs(1,:))*(X-paramEst.Xs(1,:))'/T



clf
disp('mse of Ve')
Ve = Ke*II;
VeEst = paramEst.Ke*II;


subplot(211)
plot (dt*(1:T),Ve,'r'); hold on
plot (dt*(1:T),VeAEC,'y');
plot (dt*(1:T),VeEst);
legend ('True', 'AEC', 'FIT')
xlabel('ms')
ylabel('mV')
title (sprintf('Electrode voltage V_e (MSE_{AEC} =%.2f, MSE_{FIT}=%.4f)',mseVeAEC, mseVeEMKe))
set(gca,'xlim',[0 100])
subplot(212);hold off
plot (dt*(1:T),X(1,:),'r');hold on
plot (dt*(1:T),VmAEC,'g');
plot (dt*(1:T),VmEMKe,'c--');
plot (dt*(1:T),paramEst.Xs(1,:));
legend ('True', 'AEC', 'FIT Ke','FIT smooth')
xlabel('ms')
ylabel('mV')
title (sprintf('Membrane voltage (MSE_{AEC}=%.2f, MSE_{FIT Ke}=%.2f,MSE_{FIT smooth}=%.2f)',mseVmAEC, mseVmEMKe, mseVmEMsm))
set(gca,'xlim',[0 100])

%% draw figures
subplot(612);hold off
plot (dt*(1:T), Iinj,'k--'); hold on;
plot (dt*(1:T), Ifil,'r');
plot (dt*(1:T), IfilEst);
%plot (dt*(1:T),[Iinj(:) Ifil(:)])
legend ('I_{app}', 'True I_{inj}', 'Estimated I_{inj}')
xlabel('ms')
ylabel('\muA')
title ('Input current')


%% to compare with previous results

FIT = paramEst
save est X Y II alpha beta C Ke Q R Xo Po FIT

