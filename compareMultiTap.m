
% 
% addpath ../CompareKSandEM/Kalman
% addpath ../CompareKSandEM/Kalman/jptools
addpath jptools


clear


%% basic params for estimation
ns = 1000;  % number of samples used to estimate
%ns = input ('number of samples? = ')
rep = 1;
num_trial = 8;
max_iter_per_trial = 200;


%% Choose whether to generate synthetic data or load from HH
%DATA_TYPE ='SYNTH';
%DATA_TYPE ='HH';
DATA_TYPE ='CellB';


RUN_GD = 1;


% if (input('test all Rs? (1/0)')) 
%     Rs = [0.5000    1.0000    2.0000    4.0000    8.0000]
% else	
%     Rs = 0.5
% end
Rs = 0.5

        
%% iterate for different measurement variances

for R = Rs;

    %% repeat for different input & output pairs 
    repMax = 1;
    for rep = 1:repMax

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% genrate or load data 

        switch DATA_TYPE
            case 'SYNTH'
                %% generate data on-the-fly
                dimXsynth = 3;       % number of taps = dimension of latent variable of true param
                dimY = 1;
                dimU = 1;

                %aT = rand(1,dimX); aT = aT/sum(aT)*0.9
                rr = [.95 .5 .2]';  % specify poles (roots) of IIR filter
                avec = real(poly(rr));  % coefficients of IIR filter
                aT = -avec(2:end)
                %aT = [0.7 0.15 0.1];

                beta = 0.01;
                c = 1;
                d = 1;
                Vrev = 0;
                q = 0.01;
        %         r = 0.01;
    %             r = R;
        %         q = 0.1;
        %         r = 0.25;

                Xo = zeros(dimXsynth,1)
                initPrior = 0.01;
                Po = initPrior*eye(dimXsynth,dimXsynth)

                if ~exist('ns','var')
                            ns = input('number of samples? ');
                end
                if isempty(ns)
                    ns = 5000;
                end

                Iapp = randn(1,ns);

                Ke =[0 0.8 0.4 0.1 0.1 0.05];
                U = filter(Ke, 1, Iapp);

                %         U = filter([0.25 0.25 0.25 0.25], [1 0.9], I);
                %         U = filter([0.25 0.25 0.25 0.25], 1, U);
                %         U = filter([0.25 0.25 0.25 0.25], 1, U);
                %         U = filter([0.25 0.25 0.25 0.25], 1, U);
                %         U = filter([0.25 0.25 0.25 0.25], 1, U);
                %U = sin((1:ns)/25)+ 
                %U = ones(1,ns);

                A = zeros(dimXsynth,dimXsynth);
                A(1,:) = aT;
                A(2:end,1:end-1) = eye(dimXsynth-1);
                B = zeros(dimXsynth,1);
                B(1,1) = beta;
                C = [c zeros(1,dimXsynth-1)];
                D = d;
                Q = zeros(dimXsynth,dimXsynth);
                Q(1,1) = q;
        %         R = r;

                fprintf ('generating %d samples\n', ns)
                [X,Y] = generate_lds(U, A, B, C, D, Q, R, Xo, Po);

                Iapp = Iapp';
                U = U';
                V = X(1,:)';     % true membrane voltage
                Y = Y' + Vrev;   % measured voltage

                dt = 0.1;

                %% smoothing with true param! 
                % Q. what is the minimum MSE achievable?
                XhatTrueParam = kalman_smth(Y'-Vrev, U', A, B, C, D, Q, R, Xo, initPrior);
                XhatTrueParam = XhatTrueParam(1,:);

                % plot([X(1,:)' XhatTrueParam(1,:)'])
                mseTrueParam = (X(1,2:end)-XhatTrueParam(1,2:end))*(X(1,2:end)-XhatTrueParam(1,2:end))'/(size(X,2)-1)

            case 'HH'
                %% load data    
                %dataFolderName = '../DataHH_2012-04-23_nke4_N11000_rep50'  
                dataFolderName = '.';
                filename = fullfile(dataFolderName,sprintf('HH%s_R%.2f.mat','spike',R));

                load (filename,'Iapps','Ys','Iinjs','Vms','Ke','dt')
                %% read data 

                Iapp = Iapps(1:ns,rep);     % input current 
                Y = Ys(1:ns,rep);           % measured voltage


                %% filter Iapp with Kernel Ke
                % for simplicity true Ke is used.
                % In real experiments, estimate from sub-treshold regime is used.
                KeEst = Ke; %0.4000    3.5000    1.5000    0.1500

                % input current actualy injected to the neuron
                U = filter(KeEst, 1, Iapp);

    %             Iinj = Iinjs(1:ns,rep);     % true injected current Ke*Iapp
                V = Vms(1:ns,rep);      % true membrane potential 

                Rest = R;
                betaEst = 0.014;
                q = 1;
            case 'CellB'        
                
                
                
                
                
                %% read estimated kernel, beta, and R from single tap estimates
                SUB=load('fit_CellBsub_rep1-10_N5100_itr50_trial32')
                clf;
                subplot(221)
                plot(SUB.KeEsts(:,:,1))
                title ('estimated kernel from the sub-threshold regime')
                KeEst= mean(SUB.KeEsts(:,:,1),2);
                
                subplot(222)
                betaEsts = SUB.paramEsts(2,:,1);
                betaEst = mean(betaEsts)
                errorbar(mean(betaEsts,2), std(betaEsts,[],2), 'o')
                
                title('estimated beta')
                
                subplot(223)
                Qests = SUB.paramEsts(end-1,:,1);
                errorbar(mean(Qests,2), std(Qests,[],2), 'o')
                title('estimated Q')
                set(gca,'xtick', [1 ])
                set(gca,'xticklabel',{'Q'})
                
                subplot(224)
                Rests = SUB.paramEsts(end,:,1);
                Rest = mean(Rests);
                errorbar(mean(Rests,2), std(Rests,[],2), 'o')
                title('estimated R')
                set(gca,'xtick', [1])
                set(gca,'xticklabel',{'R'})
                
                set(gcf, 'paperposition', [0 0 8 6]) 
                set(gcf, 'papersize', [8 6])
                saveas(1,'estimation_sub-threshold.pdf')

                %% read spiking data 
                filename = 'Cell_04_15_2010_BD_spike.mat'

                load (filename,'Iapps','Ys','Y2s','dt')
                
                Iapp = Iapps{3}(1:ns);     % input current 
                Y = Ys{3}(1:ns);  
                Y2 = Y2s{3}(1:ns);  
                
                %% input current actualy injected to the neuron
                U = filter(KeEst, 1, Iapp);
                
                
        end

        % %     %% save partial data (for generating smaller data set)
        % %     Iapps = Iapps(:,1:10);
        % %     Ys = Ys(:,1:10);
        % %     Iinjs = Iinjs(:,1:10);
        % %     Vms = Vms(:,1:10);
        % %     
        % %     filename = sprintf('HH%s_R%.2f.mat','spike',R)
        % %     save (filename,'Iapps','Ys','Iinjs','Vms','Ke','dt')

    
        %% repeat with different number of taps 
        %for dimX = 1:5       % number of taps = dimension of latent variable of true param
        for dimX = 2       % number of taps = dimension of latent variable of true param
       
    %         disp(sprintf('Fit %dth course',rep))

            %% find params from linear regression => used for init of EM & GD
            disp('initial guess from linear regression')
            MM = fliplr(makeStimRows([0;Y(1:end-1)-U(1:end-1)],dimX)); % make design matrix for regressing Y on its history
            MM = [MM ones(size(MM,1),1)];  % mean is also estimated!
            aT2Vo2=MM\(Y-U);             % solve regression as initial values for alphas
            aT2 = aT2Vo2(1:end-1);        % use residual variance to initialize estimate for Q
            disp('corresponding roots are')
            disp(roots([1;-aT2]))
            Vo2 = aT2Vo2(end);
            residue = var(Y-U-MM*aT2Vo2);
            Q2 = residue-R;  % use residual variance to initialize estimate for Q 
            if Q2>10
                fprintf('[WARNING] initial Q=%.1f is too large. Replaced it with 10.', Q2);
                Q2 = min(Q2,10);
            end
            Q2 = max(Q2,1e-9);
            


            % initial value is selected uniformly from the given interval
            % range for the vector alpha (At)
            AtLow = aT2-abs(aT2*0.1)-eps;
            AtHigh = aT2+abs(aT2*0.1)+eps;
            % range for [beta V_rev Q R]
            % use true beta and R
            ParamLow =[betaEst, Vo2+5, Q2/2, Rest];
            ParamHigh =[betaEst, Vo2-5, Q2*2, Rest];
    % %         % use noisy R  (Test)
    % %         ParamLow =[1e-3, mean(Y)-mean(U)+5, 1, R*0.99];
    % %         ParamHigh =[5e-2, mean(Y)-mean(U)-5, 10, R*1.01];


            initPrior = 10;     % prior: variance of V at time 0


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% 1. Call EM routine, which does the following
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 1) initialize params and give back to initAt and initParam
            % 2) repeat estimation num_trial times and store estimates in all
            % trials into allAt and alllParam
            % 3) return final estimate (with largest LL) in EM.AtEst and paramEstMulti
            fprintf ('estimating by EM\n')
            [EM.AtEst, EM.paramEstMulti, EM.LL, EM.LLs, initAt, initParam, allAt, allParam] = estimParamEM_multi(Y, U, AtLow, AtHigh, ParamLow, ParamHigh, initPrior, max_iter_per_trial, num_trial);
            LLtrials = max(EM.LLs);
            itrTrials = sum(~isnan(EM.LLs))

            %% Smoothe with parameters estimated by EM
            A = zeros(dimX,dimX);
            A(1,:) = EM.AtEst;
            A(2:end,1:end-1) = eye(dimX-1);
        %                    B = zeros(dimX,1);
        %                     B(1,1) = paramEstMulti(1);      % beta
            %B = EM.paramEstMulti(1);      % beta
            B = [EM.paramEstMulti(1); zeros(dimX-1,1)];      % beta
            C = [1 zeros(1,dimX-1)];
            D = 1;
            Q = zeros(dimX,dimX);
            Q(1,1) = EM.paramEstMulti(3);
            Rhat = EM.paramEstMulti(4);
            Xo = zeros(dimX,1);

            disp('smoothing with multi-tap results')
            if (size(Y,1)==1)
                [Xs Ps Pcs] = kalman_smth(Y-EM.paramEstMulti(2), U, A, B, C, D, Q, Rhat, Xo, initPrior);
            else 
                [Xs Ps Pcs] = kalman_smth(Y'-EM.paramEstMulti(2), U', A, B, C, D, Q, Rhat, Xo, initPrior);
            end
            EM.Vest = Xs(1,:)' + EM.paramEstMulti(2);

            if exist('V','var')
                EM.mseMulti(rep) = (V(2:end)'-EM.Vest(2:end)')*(V(2:end)-EM.Vest(2:end))/(length(V)-1)
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% plot EM results
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % find  peak and plot around it 
            [~, idxMax] = max(EM.Vest);
            idx = max(idxMax-100/dt,1):min(idxMax+100/dt,ns);
            clf;
            numRow=6;

            % plot trajectories 
            subplot(numRow,1,1);
            if exist('Iapp','var')
                plot(dt*(idx),Iapp(idx),'r');hold on
    %             plot(dt*(idx),Iinj(idx));
                legend('I_{app}','I_{inj}'); legend boxoff
                xlabel('ms')
                ylabel('uA/cm^2')
                title ('applied current')
                box off
            end

            subplot (numRow,1,2)
            plot (dt*(idx),Y(idx))
            if exist('Y2s','var')    
                hold on;
                plot (dt*(idx),Y2(idx),'g')
                legend ('measure', 'reference');
            end
            title ('measured voltage (V_{rec})')
            box off
            ylabel('mV')


            subplot(numRow,1,3);
            plot(dt*(idx),EM.Vest(idx)); hold on 
    %         plot(dt*(idx),V(idx),'r--');
    %         legend('EM','true'); legend boxoff
            if exist('EM.mseMulti','var')
                title(sprintf('Membrane voltage (V_m) MSE(EM)=%.4f',EM.mseMulti(rep)));
            end


            xlabel('ms')
            ylabel('mV')
            box off



            subplot(numRow,2,7)
            plot(EM.LLs)
            set(gca,'xlim', [1 max(itrTrials)])
            title ('LL(EM)')
            xlabel ('iteration')
            box off

            subplot(numRow,2,8)
            plot(2:(size(EM.LLs,1)),diff(EM.LLs))
            set(gca,'xlim', [1 max(itrTrials)])
            title ('diff(LL(EM))')
            xlabel ('iteration')
            box off


            subplot(numRow,2,9)
            %hist(LLtrials)
            plot (LLtrials, 'o--')

            % high light max LL
            hold on;
            [LLmax maxTrial] = max(LLtrials);
            plot(maxTrial,LLmax,'ro')

            title ('final LLs (EM)')
            xlabel ('trial')
            box off

            subplot(numRow,2,10)
            plot (itrTrials, 'o--'); hold on
            plot(maxTrial,itrTrials(maxTrial),'ro')
            title ('iterations before convergence')
            ylim = get(gca,'ylim');
            set(gca,'ylim',[0 ylim(2)])
            %hist(itrTrials)
            %title ('histogram of iteration before convergence')
            xlabel ('trial')
            box off



            
            

            subplot(numRow,4,21);
            plot(initAt','g+:'); hold on
            plot(allAt','bo:'); 
            plot(maxTrial, allAt(:,maxTrial)','ro'); 
            xlabel('trial')
            title ('aT')
            box off

            subplot(numRow,4,22);
            root = roots([1; -EM.AtEst(:)]); hold on
            plot(real(root),imag(root), 'or')
            title ('roots of aT')


            subplot(numRow,4,23);
            plot(initParam(2,:)','g+:'); hold on
            plot(allParam(2,:)','bo:'); hold on
            plot(maxTrial, allParam(2,maxTrial)','ro'); hold on
            xlabel('trial')
            title ('V_{rev}')
            box off

            subplot(numRow,4,24);
            plot(initParam(3,:)','g+:'); hold on
            plot(allParam(3,:)','bo:'); hold on
            plot(maxTrial, allParam(3,maxTrial)','ro'); hold on
            xlabel('trial')
            title ('Q')
            set(gca,'yscale','log')
            box off

            legend ('initial','EM', 'EM with max LL'); legend boxoff


            drawnow

%% 
if RUN_GD
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Now, call gradient descent based method by JP
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            switch DATA_TYPE
                case 'SYNTH'
                    %prstrue = [aT, log(1/q), log(1/R), Vrev, log(beta)]';
                    prstrue = [aT; log(1/q); Vrev];
                case 'HH'
            end

            fprintf ('estimating by GD\n')
            prshat = zeros (dimX+2,num_trial);
            GD.LLs = nan*ones(num_trial,1);
            for trial = 1:num_trial
                % start from the same init as EM
                %prs0 = prstrue+randn(size(prstrue))*.1;

                % Alter initial values so roots are all abs<1
                aT0 = initAt(:,trial);
                rts = roots([1;-aT0]);
                if any(abs(rts)>=1)
                    rts = rts./((max(abs(rts))+.01));
                end
                aT1 = poly(rts);
                aT1 = -aT1(2:end)';

                % toggle between ways of setting initial params for alphas and Q

                if 0  % use the vals from EM, rounded so alpha roots are <1
                    prsInit = [aT1; log(1./initParam(3,trial))];
                else
                    M = fliplr(makeStimRows([0;Y(1:end-1)],dimX)); % make design matrix for regressing Y on its history
                    aT2 = M\Y;  % solve regression as initial values for alphas
                    Q2 = var(Y-M*aT2);  % use residual variance to initialize estimate for Q
                    prsInit = [aT2; log(1./Q2)];
                end
                % embed these params in the parameter vector
                prs0 =  [prsInit; initParam(2,trial)];  %[aT -logQ mu_y]

                % define loss function
                lFunc = @(prs)(neglogmargLi_LDSvec_DynamPrsOnly(prs, Y(:), dimX,U(:), Rest, betaEst));   % updated to fix beta and r
                
                opts = optimset('display', 'iter','largescale','off','maxiter',2000,'maxfuneval',25000,'tolfun',1e-16,'tolx',1e-16);
                [prs negLL] = fminunc(lFunc,prs0, opts); % compute estimate
                
                % store estimated param
                prshat(:,trial) = prs;

                % store LL
                GD.LLs(trial) = -negLL;
                
            end

            %prshat = repmat(prstrue,1,num_trial);


            %% choose best estimate (with largest LL) and smoothe
            % replace -1e25 with NAN
            idxNan = GD.LLs < -1e24;
            GD.LLs(idxNan) = NaN;




            [mm maxTrialGD] = max(GD.LLs)
            GD.AtEst = prshat(1:dimX,maxTrialGD)

            %% Smoothe with parameters estimated by EM
            A = zeros(dimX,dimX);
            A(1,:) = prshat(1:dimX,maxTrialGD)
            A(2:end,1:end-1) = eye(dimX-1);
        %                    B = zeros(dimX,1);
        %                     B(1,1) = paramEstMulti(1);      % beta
            %B = EM.paramEstMulti(1);      % beta
            B = [betaEst; zeros(dimX-1,1)];      % beta
            C = [1 zeros(1,dimX-1)];
            D = 1;
            Q = zeros(dimX,dimX);
            Q(1,1) = exp(-prshat(dimX+1,maxTrialGD));
            Rhat = R;
            Xo = zeros(dimX,1);

            disp('smoothing with multi-tap result of GD ')
            if (size(Y,1)==1)
                [Xs Ps Pcs] = kalman_smth(Y-prshat(dimX+2,maxTrialGD), U, A, B, C, D, Q, Rhat, Xo, initPrior);
            else 
                [Xs Ps Pcs] = kalman_smth(Y'-prshat(dimX+2,maxTrialGD), U', A, B, C, D, Q, Rhat, Xo, initPrior);
            end
            GD.Vest = Xs(1,:)' + prshat(dimX+2,maxTrialGD);

            if exist('V','var')
                GD.mseMulti(rep) = (V(2:end)'-GD.Vest(2:end)')*(V(2:end)-GD.Vest(2:end))/(length(V)-1);
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Plot JP results 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     ahat = prshat(1:dimX,:);
        %     sigxhat = sqrt(exp(-prshat(dimX+1)));
        %     sigyhat = sqrt(exp(-prshat(dimX+2)));
        %     muyhat = prshat(dimX+3);
        %     betahat = prshat(dimX+4:end);


            LINE_TYPE_JP = 'sk--'
            LINE_TYPE_JP_MAX_LL = 'sr'

            % plot estimated membrane potential from GD with true 
            subplot(numRow,1,3);
            plot(dt*(idx),GD.Vest(idx),'k-.'); hold on 
            if exist('V','var')
                plot(dt*(idx),V(idx),'r--');
                legend ('EM','GD','true')
            end
            if exist('Y2','var')
                plot(dt*(idx),Y2(idx),'g'); hold on 
                legend ('EM','GD','referene')
            end
            
            if exist('XhatTrueParam','var')
                plot(dt*(idx),XhatTrueParam(idx),'r--');
                title(sprintf('Membrane voltage (V_m) MSE(EM)=%.4f, MSE(GD)=%.4f, MSE(true)=%.4f',EM.mseMulti(rep),GD.mseMulti(rep),mseTrueParam));
            else
                switch DATA_TYPE
                case {'SYNTH', 'HH'}
                    title(sprintf('Membrane voltage (V_m) MSE(EM)=%.4f, MSE(GD)=%.4f',EM.mseMulti(rep),GD.mseMulti(rep)));
                end
                
            end
            


            % plot final LLs for each trial
            subplot(numRow,2,9)
            plot (GD.LLs, LINE_TYPE_JP)
            plot (maxTrialGD,GD.LLs(maxTrialGD), LINE_TYPE_JP_MAX_LL)

            legend ('EM','EM max','GD','GD max'); legend boxoff


            subplot(numRow,4,21);   % aT
            plot(prshat(1:dimX,:)',LINE_TYPE_JP)
            plot(maxTrialGD,prshat(1:dimX,maxTrialGD)',LINE_TYPE_JP_MAX_LL)

            subplot(numRow,4,22);
            root = roots([1; -prshat(1:dimX,maxTrialGD)]);
            plot(real(root),imag(root), LINE_TYPE_JP_MAX_LL)
            title ('roots of aT')
            
            subplot(numRow,4,23);   % Vrev
            plot(prshat(dimX+2,:),LINE_TYPE_JP)
            plot(maxTrialGD,prshat(dimX+2,maxTrialGD)',LINE_TYPE_JP_MAX_LL)
            subplot(numRow,4,24);   % Q
            plot(exp(-prshat(dimX+1,:)),LINE_TYPE_JP)
            plot(maxTrialGD,exp(-prshat(dimX+1,maxTrialGD)),LINE_TYPE_JP_MAX_LL)
            
end
            %% Plot true param for synthetic data case 
            switch DATA_TYPE
                case 'SYNTH'
                            
                    LINE_TYPE_TRUE = 'k--x'
                    subplot(numRow,5,26);   % aT
                    plot(repmat(aT,num_trial),LINE_TYPE_TRUE); hold on        
                    subplot(numRow,5,27);   % beta
                    plot(repmat(beta,num_trial),LINE_TYPE_TRUE); hold on        
                    subplot(numRow,5,28);   % beta
                    plot(repmat(Vrev,num_trial),LINE_TYPE_TRUE); hold on        
                    subplot(numRow,5,29);   % q
                    plot(repmat(q,num_trial),LINE_TYPE_TRUE); hold on        
                    subplot(numRow,5,30);   % beta
                    plot(repmat(R,num_trial),LINE_TYPE_TRUE); hold on        

                    legend ('initial','EM','EM-max','GD','GD-max','true'); legend boxoff
                    
                    set(gcf, 'paperposition', [0 0 12, 18]) 
                    set(gcf, 'papersize', [12 18]) 
                    filename = sprintf('%s_EM_vs_GD_dim%d_Q%.2f_R%.2f_rep%d.pdf',DATA_TYPE,dimX,q,R,rep);   % we know q for synthetic case
                    saveas(1,filename)


                    filename = sprintf('%s_EM_vs_GD_dim%d_Q%.2f_R%.2f_rep%d.mat',DATA_TYPE,dimX,q,R,rep)
                    save(filename)
            
                case 'HH'
                    legend ('initial','EM','EM-max','GD','GD-max'); legend boxoff
                    
                    set(gcf, 'paperposition', [0 0 12, 18]) 
                    set(gcf, 'papersize', [12 18]) 
                    filename = sprintf('%s_EM_vs_GD_dim%d_R%.2f_rep%d.pdf',DATA_TYPE,dimX,R,rep)
                    saveas(1,filename)


                    filename = sprintf('%s_EM_vs_GD_dim%d_R%.2f_rep%d.mat',DATA_TYPE,dimX,R,rep)
                    save(filename)
            end





            

        end
    end  % end for dimX
end   % end of for Rs



return 


%% Root analysis 
numRow = 5;
numCol = length(Rs);

% unit circle;
theta = 0:.01:(2*pi);
UC = [cos(theta); sin(theta)];

clf;
mseEM = zeros(5,length(Rs));
qHatEM = zeros(5,length(Rs));
for dim = 1:5
    for i = 1:length(Rs);
        R=Rs(i);
        %filename = sprintf('%s_EM_vs_GD_dim%d_Q%.2f_R%.2f_rep%d.mat',DATA_TYPE,dim,q,R,rep)
        filename = sprintf('%s_EM_vs_GD_dim%d_R%.2f_rep%d.mat',DATA_TYPE,dim,R,rep)
        load(filename, 'EM', 'GD', 'aT','mseMulti')
        
        % calc and plot root
        rootEM = roots([1 -EM.AtEst]);
        if exist('GD','var')
            rootGD = roots([1; -GD.AtEst]);
        end
        if exist('aT','var')
            rootTrue=roots([1 -aT]);
        end
        

        subplot(numRow,numCol,i+(dim-1)*numCol)
        plot(real(rootEM),imag(rootEM),'ob'); hold on
        if exist('GD','var')
            plot(real(rootGD),imag(rootGD),'sg');
        end
        if exist('aT','var')
            plot(real(rootTrue),imag(rootTrue),'xr');
        end
        
        % plot unit circle
        plot (UC(1,:), UC(2,:), ':k')
        box off
        
        % MSE
        mseEM(dim,i) = EM.mseMulti;
        
        % check q
        qHatEM(dim,i) = EM.paramEstMulti(end-1);
        
        if i==1 & dim==1
            if exist('GD','var')
                if exist('aT','var')
                    legend ('EM','GD','true'); %legend boxoff
                else
                    legend ('EM','GD'); %legend boxoff
                end
                
            else 
                if exist('aT','var')
                    legend ('EM', 'true')
                end
            end
        end
        
        if dim==1
            title (sprintf('R=%.1f',R))
        end
        
        if i==1
            ylabel (sprintf('tap=%d',dim))
        end

    end
end

set(gcf, 'paperposition', [0 0 15, 12]) 
set(gcf, 'papersize', [15 12]) 
filename = sprintf('%s_EM_roots.pdf',DATA_TYPE)
saveas(1,filename)


%%
clf
subplot(121)
plot(log2(Rs),qHatEM','o-'); xlabel('log2(R)')
if exist('q','var')
    hold on
%    plot(log2(Rs),repmat(q,1,length(Rs)),'k--+')
end
legend('tap=1','tap=2','tap=3','tap=4','tap=5','true','location','nw');legend boxoff
title ('Qhat'); box off
subplot(122)
plot(log2(Rs),log2(mseEM(1:3,:)'),'o-'); xlabel('log2(R)')
hold on; plot(log2([Rs(1) Rs(end)]),log2([Rs(1) Rs(end)]),'k+--')
legend('tap=1','tap=2','tap=3','location','nw');legend boxoff
%axis equal
title ('log2(MSE)'); box off



set(gcf, 'paperposition', [0 0 8, 3]) 
set(gcf, 'papersize', [8 3]) 
filename = sprintf('%s_EM_q_mse.pdf',DATA_TYPE)
saveas(1,filename)
                    
%%
clf;
subplot(211)
imagesc(1:5,Rs,qHatEM); xlabel('R');ylabel('tap');colormap gray;colorbar
title ('Qhat')
subplot(212)
imagesc(1:5,Rs,mseEM); xlabel('R');ylabel('tap');colormap gray;colorbar
title ('MSE')