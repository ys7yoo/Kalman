%% code copied from src.vis/EMKalman_0606/fit_EM/fit_EM_Rs.m
% which has the latest code for HH

clear
close all

addpath('../AEC')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% chose data & method to fit 
%EXP.DATA = 'LDS';
%EXP.DATA = 'LDS_S';  % small scale 
%EXP.DATA = 'HH';
%EXP.DATA = 'CellB'
% EXP.DATA = 'Cell2'
switch input('Data to fit (0:LDS, 1:HH, 2:CellB, 3:Cell2)? ')
    case 1
        EXP.DATA = 'HH'
    case 2
        EXP.DATA = 'CellB'
    otherwise
        return
end

EXP.MODE = input('MODE (0:SUB,1:SPIKE-SINGLE,2:SPIKE-MULTI,3:AEC)? ')
switch EXP.MODE
    case 0
        EXP.METHOD = 'sub';
    case 1
        EXP.METHOD = 'spike-single'    %spike data but try to fit with single tap.
    case 2
        EXP.METHOD = 'spike'        
    case 3
        EXP.METHOD = 'AEC'        
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for EM algorithm (Single tap)
ESTIMATE_SINGLE_TAP_EM=1;

EM.INIT_KE_BY_AEC=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for EM algorithm (multi tap)
ESTIMATE_MULTI_TAP_EM=1     % use the same param to smoothe spiking data
%ESTIMATE_MULTI_TAP_EM=0    

addpath('../Kalman')
addpath('../Kalman/jptools')
addpath('../YsTools')


if EXP.MODE~=3
    EM.num_trial = 32;
    
    %EM.max_iter = 50;
    EM.max_iter = 100;
    %EM.max_iter = 500; 
    %EM.max_iter = 5000; 
    %EM.max_iter = 1e6;  % Huge. Converges for all Rs BUT, not practical params
    %EM.max_iter = input ('Max iteration per trial? ');



    % choose M-step method
    EM.MstepMethod = 0;   % partial derivative 
    %EM.MstepMethod = 1;   % matlab opt function (fminunc)
    %EM.MstepMethod = 2;   % matlab opt with constraint (fmincon)
    %EM.MstepMethod = -1;   % compare method (for debuggin)

    % to add constraint or not
    EM.MstepConstraint = 0;
    if EM.MstepMethod == 2
        EM.MstepConstraint = 1;
    end

    EM.checkConverged = 1;
    EM.checkDecreased = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% params for gradient descent method 
ESTIMATE_MULTI_TAP_GD=1
GD.num_trial = EM.num_trial;
%GD.num_trial = 2; % for debug
addpath ../Kalman/jptools/




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load some basic params depending on data set 
switch EXP.DATA
    case 'LDS'
%         Qs = .1;
%         Rs = .1;

%         Qs = 1
%         Rs = 1

        Qs = .1
%         Qs = 1
        Rs = 10.^(-1:.25:0)        % for batch mode
        
        filename = sprintf('../Data/LDS_Q%.2f_R%.2f.mat',Qs(1),Rs(1));
        load (filename, 'numRep', 'param');
        paramTrue = param;
        EXP.numRep = numRep;
        EXP.numRep = 100;
        %EXP.numRep = 1;
        EM.M = length(paramTrue.Ke);
        
        
        EM.T=5000;
        EXP.TRUE = 1;
        dt = .1;
    case 'LDS_S'
        filename = sprintf('../Data/LDS_Q%.2f_R%.2f.mat',Qs(1),Rs(1));
        load (filename, 'numRep', 'param');
        paramTrue = param;
        EXP.numRep = 1;
        EM.M = length(paramTrue.Ke);
        

        EM.T=1000;
        EM.num_trial = 1;
        EXP.TRUE = 1;
        dt = .1;        
        
    case 'HH'   % fit HH model
        Qs = 1
        %Rs = 10.^(-1:.25:0)        % for batch mode
        %Rs = 10.^(-2:.25:-1)        % for batch mode
        %Rs = 10.^(-2:.25:0)         % for now, stick to this range (Mar. 22, 2012)
        % Rs = 10.^(-1:.5:1)        % readjust noise range larger noise is going to be more interesting
        Rs=2.^(-1:3)
        
%         dataFolderName = '../DataHH_20120315_ext0';
%         dataFolderName = '../DataHH_20120315_ext1000';
%         dataFolderName = '../DataHH_2012-03-15_nke8_N10000_ext0'
%         dataFolderName = '../DataHH_2012-03-15_nke8_N10000_ext1000';


        % after moving subthreshold data (03-16) and spiking data (04-23)
        % into a folder
        
        % less spiking data Idc = 0 Apr 23, 2012
%         dataFolderName = '../DataHH_2012-04-23_nke4_N11000_rep50';
        dataFolderName = '../DataHH_2012-05-02_nke4_N21000';        % long data on melville
        
        if exist(dataFolderName,'dir')
            numRepMax = input ('numRepMax = ? ');
        else
            sprintf('folder %s doesn''t exist. Try to load smaller set', dataFolderName)
        

            switch EXP.METHOD
                case {'sub'}
                    % go back to shorter ke
                    %dataFolderName = '../DataHH_2012-03-16_nke5_N5100_ext0'   % has
                    %zero in the first tap
                    %dataFolderName = '../DataHH_2012-03-16_nke4_N5100_ext0'

                    dataFolderName = '../DataHH_2012-03-16_nke4_N1100'      % small data for debugging       
                    numRepMax = 1;

                case {'spike','spike-single'}
                    % less spiking data Idc = 0 Apr 23, 2012
                    dataFolderName = '../DataHH_2012-04-23_nke4_N11000_rep50'       
                    numRepMax = 1;  


            end
		
%         Rs = 10.^(-1:1)        % even coarser sampling
        end


        % data to load
        R = Rs(1);
        switch EXP.METHOD
            case {'sub'}
                filename = fullfile(dataFolderName,sprintf('HH%s_R%.2f.mat','sub',R));
        
            case {'spike','spike-single'}
                filename = fullfile(dataFolderName,sprintf('HH%s_R%.2f.mat','spike',R));
                
            otherwise  % general version for the future
                filename = fullfile(dataFolderName,sprintf('HH%s_R%.2f.mat',EXP.METHOD,R));
        end

        load (filename, 'numRep', 'Ke');
        %EXP.numRep = numRep;
        %EXP.numRep = size(Ys,2);
        if exist('numRep','var')
            EXP.numRep = min(numRepMax,numRep);
        else
            EXP.numRep = 1;
        end
        EM.M = length(Ke);
        
        
        EM.T=5000;
        %EM.T = input('Number of samples to use for estimation? ');
        
        
        
        EXP.TRUE = 1;
        dt = .1;
     
        
    case 'CellB'
        EXP.numRep = 10;
        %EXP.numRep = 1;

        dt = 0.1;
        EM.M = 5/dt;
        %EM.M = 25;
        
        
        EXP.TRUE = 0;
        Qs = 1;
        
        Rs = [0.1 0.25 0.5 0.75];   
        % real data R is the variance of input in the linear regime.
        % for multi-tap fitting, R determins the linear estimate to use 
        
	
	
        %EM.T = 100;
     	EM.T = 5000;
    % 	EM.T = 10000;
        %EM.T=length(Ys); fit full data
        
        %dt = 80/1000;
    case 'Cell2'
        EXP.numRep = 10;
        %EM.M = 50;
        EM.M = 25;
        
        
        EXP.TRUE = 0;
        Qs = 1;
        Rs = 1;
	
	
        %EM.T = 100;
     	EM.T = 5000;
    % 	EM.T = 10000;
        %EM.T=length(Ys); fit full data
        
        %dt = 80/1000;        
end


%EM.Next = input('Number of initial samples to exclude for evaluation? ')
EM.Next = 100;
EM.T = EM.T + EM.Next;

        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory for results
% and read single-tap results for multi-tap estimation
numR = length(Rs);
mse = NaN*ones (EXP.numRep, numR);
Vests = NaN*zeros(EM.T,EXP.numRep,numR);

LLfinals = NaN*ones (EXP.numRep, numR);
maxItrs = NaN*ones (EXP.numRep, numR);
switch EXP.METHOD
    case {'sub','spike-single'}
        
        paramEsts = NaN*zeros (5,EXP.numRep, numR);
        KeEsts = NaN*zeros (EM.M,EXP.numRep, numR);
        
        mseByAEC = NaN*zeros (EXP.numRep, numR);
        mseEMKeLin = NaN*zeros (EXP.numRep, numR);
        mseKe = NaN*zeros (EXP.numRep, numR);

        
        VestsByAEC = NaN*zeros(EM.T,EXP.numRep,numR);
        
        dimX = 1;   % single tap means dimension of hidden state is 1
        
    case {'spike','post-lin'}
        
        dimX = input('dimX=? ')

        %% initialize memory for estimate resuls 
        
        % linear regressions with multitap from true V
        paramEstTrueV = NaN*zeros (dimX+4,EXP.numRep, numR);
        mseTrueV = NaN*zeros (EXP.numRep, numR); %mse for smoothing with regression wit true V
        
        % linear fit with multitap from the measurement 
        paramEstLinRegMulti = NaN*zeros (dimX+2,EXP.numRep, numR);   % store aT, Vo, and Q
        mseLinRegMulti = NaN*zeros (EXP.numRep, numR);
        
        % smoothing with single-tap results
        mseSingleTap = NaN*zeros (EXP.numRep, numR);
        
        % re-fit with EM
        if (ESTIMATE_MULTI_TAP_EM)
            AtEsts = NaN*zeros (dimX,EXP.numRep, numR);
            paramEstMultis = NaN*zeros (4,EXP.numRep, numR);
            mseMultiTap = NaN*zeros (EXP.numRep, numR);
            
        end
        
        %% load params from single tap fitting
        switch EXP.DATA
            case {'LDS','LDS_S'}
                %??
            
            
            case 'HH'
        
               foldernameSingleTap = '.';   % generally load single-tap result in the same folder
               %foldernameSingleTap = '../2012-05-17_HH_nke4_additional'; % or fix manually
               if exist(foldernameSingleTap,'dir')   % fix manually 
                   %filenameSingleTap = fullfile(foldernameSingleTap, sprintf('fit_HHsub_rep1-%d_N%d_itr500_trial32.mat',EXP.numRep,EM.T))
                   %filenameSingleTap = fullfile(foldernameSingleTap, sprintf('fit_HHsub_rep1-50_N%d_itr500_trial32.mat',EM.T))
                   %filenameSingleTap = fullfile(foldernameSingleTap, 'fit_HHsub_rep1-50_N5100_itr500_trial32.mat'); % fix N as well
                   filenameSingleTap = fullfile(foldernameSingleTap, sprintf('fit_HHsub_rep1-%d_N5100_itr%d_trial32.mat',EXP.numRep,EM.max_iter)); % fix N as well

        %             % load KeEsts from here 
        %             filenameSingleTap = sprintf('../2012-05-01_HH_nke4_single_tap_detail/fit_HHsub_rep1-10_N%d_itr1000000_trial32.mat',EM.T)
        %             filenameSingleTap = sprintf('../2012-05-01_HH_nke4_single_tap_detail/fit_HHsub_rep1-10_N%d_itr2000_trial32.mat',EM.T)

        %       WHYT DID I USED THE FOLLOWING?
        %         if exist('../2012-04-25_HH_nke4_spike_single_tap2','dir')     % load  from linear fitting 
        %             filenameSingleTap = '../2012-04-25_HH_nke4_spike_single-tap2/fit_HHspike-single_rep1-50_N5100_itr5000_trial32.mat'

                else % ask 
                    % choose filename from input 
                    maxIterSingleTap = input ('max iteration of single tap? ');
                    numRepSingleTap = input ('num rep of single tap?');     % numRep of single tap result could be larger
                    %filenameSingleTap = sprintf('fit_HHsub_rep1-%d_dim1_N%d_itr%d_trial32.mat',EXP.numRep,EM.T,maxIterSingleTap)
                    %filenameSingleTap = sprintf('fit_HHsub_rep1-%d_N%d_itr%d_trial32.mat',EXP.numRep,EM.T,maxIterSingleTap)
                    filenameSingleTap = sprintf('fit_HHsub_rep1-%d_N%d_itr%d_trial32.mat',numRepSingleTap,EM.T,maxIterSingleTap)            

               end
               
            case 'CellB'
                %fit_CellBsub_rep1-10_N5100_itr50_trial32.mat
                numRepSingleTap = EXP.numRep;
                maxIterSingleTap = EM.max_iter;
                filenameSingleTap = sprintf('fit_%ssub_rep1-%d_N%d_itr%d_trial32.mat',EXP.DATA,numRepSingleTap,EM.T,maxIterSingleTap)            
            
        end
        
        
       % load here 
        singleTap = load (filenameSingleTap,'paramEsts','KeEsts','LLfinals')
            
        
        
    
    case 'AEC'  % I want to AEC and store separately
        mseKe = NaN*zeros (EXP.numRep, numR);
        
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% main iteration starts from here 
disp(EM)
disp(GD)




timeDuration = zeros(numR,1)
for i=1:numR;
%for i=6:9;
    tic;
    Q = Qs;
    R = Rs(i);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Load data 
    switch EXP.DATA
        case {'LDS','LDS_S'}
            % load data generated by HH
            filename = sprintf('../Data/LDS_Q%.2f_R%.2f.mat',Q,R);
            load (filename, 'Iapps', 'Iinjs','Vms', 'Ves', 'Ys', 'param')
            
            
        case 'HH'
            % load data generated by HH   
            disp(sprintf('Fit time courses for R=%.2f,',R))

            
            switch EXP.METHOD
                case {'sub', 'AEC'}
                    EXP.baseFilename = sprintf('HH%s_R%.2f.mat','sub',R);

                case {'spike','spike-single'}
                    EXP.baseFilename = sprintf('HH%s_R%.2f.mat','spike',R)
                    
                otherwise
                    error ('input baseFilename is not set')
%                 otherwise  % general version for the future
%                     filename = fullfile(dataFolderName,sprintf('HH%s_R%.2f.mat',EXP.METHOD,R))
            end
            
            EXP.filename = fullfile(dataFolderName,EXP.baseFilename)
            
            %filename = fullfile(dataFolderName,sprintf('HH%s_R%.2f.mat','sub',R));     
            load (EXP.filename, 'Iapps', 'Iinjs','Vms', 'Ves', 'Ys', 'tt', 'Ke')
            
            paramTrue.Cm = 1.0;
            paramTruel.gL = 0.3;
            paramTrue.Vo = 10.6;
            paramTrue.R = R;
        
            Ke
            
        case 'CellB'
            % filename to load 

            switch EXP.METHOD
                case {'sub','AEC'}
                    EXP.baseFilename = sprintf('Cell_04_15_2010_BD_n%.2f',R);
                case {'spike','spike-single'}
                    EXP.baseFilename = sprintf('Cell_04_15_2010_BD_spike');
            end
            EXP.filename = sprintf('../DataCell1/%s.mat',EXP.baseFilename)
            
                    
            % old version, manually fix 
            %filename = sprintf('../Data/cell_4_15_2010_B_80_run0p1.mat') %=> bad
            %filename = sprintf('../Data/cell_4_15_2010_B_80_run0p5.mat')
            %filename = sprintf('../Data/cell_4_15_2010_B_80_run0p25.mat')
            %filename = sprintf('../Data/cell_4_15_2010_B_80_run0p25_s2.mat')
            %filename = sprintf('../Data/cell_4_15_2010_B_80_run0p75.mat')
            %load (filename);
            
            
            

% %             DATA_SET = input('data set? = ')
% %             switch DATA_SET
% %                 % old datat set with various inject current levels.
% %                 case 1
% %                     EXP.baseFilename = 'Cell_04_15_2010_BD_n0.1';
% %                     EXP.filename = sprintf('../DataCell1/%s.mat',EXP.baseFilename)
% %                 case 2
% %                     EXP.baseFilename = 'Cell_04_15_2010_BD_n0.25';
% %                     EXP.filename = sprintf('../DataCell1/%s.mat',EXP.baseFilename)
% %                 case 3
% %                     EXP.baseFilename = 'Cell_04_15_2010_BD_n0.5';
% %                     EXP.filename = sprintf('../DataCell1/%s.mat',EXP.baseFilename)
% %                 case 4
% %                     EXP.baseFilename = 'Cell_04_15_2010_BD_n0.75';
% %                     EXP.filename = sprintf('../DataCell1/%s.mat',EXP.baseFilename)
% %                 case 5
% %                     % % 'Cell2_uncomp'
% %                     EXP.baseFilename = 'AD_uncomp_n0.1';
% %                     EXP.filename = sprintf('../DataCell2/%s.mat',EXP.baseFilename)
% %             end
            load (EXP.filename);

            if iscell(Iapps)
                Iapps = cell2mat(Iapps);
                Y2s = cell2mat(Y2s);
                Ys = cell2mat(Ys);
            end

            
% %             %% show linear kernel
% %             II=stackCols(Iapps(:,1),EM.M,0);
% %             II = [II; ones(1,size(II,2))];
% %             KlinVoLin = II'\Ys(:,1);
% %             Klin = KlinVoLin(1:end-1);
% %             VoLin=KlinVoLin(end);
% %             clf;
% %             plot(Klin)
% %             title ('Composite kernel')
% %             saveas(1, sprintf('%s_composite_kernel_M%d.pdf',EXP.baseFilename,EM.M))
% %             
% %             if 1==0
% %                 %% compare different number of taps
% %                 clf;
% %                 subplot(121)
% %                 plot(Klin50,'r'); hold on
% %                 plot(Klin500,'g--');
% %                 plot(Klin1000,'b-.');
% %                 legend ('50 lag', '500 lag')
% %                 title ('Composite kernel')
% % 
% %                 subplot(122)
% %                 semilogx(Klin50,'r'); hold on
% %                 semilogx(Klin500,'g--');
% %                 semilogx(Klin1000,'b-.');
% % 
% %                 saveas(1, sprintf('Composite kernels.pdf'))
% %             end
            
            %%
            paramTrue = [];
        case 'Cell2'
           
            filename = sprintf('../Data2/AD_uncomp_n0.1.mat');  % (almost) linear data
            load (filename);
            

            
% %             %% show linear kernel
% %             II=stackCols(Iapps(:,1),EM.M,0);
% %             II = [II; ones(1,size(II,2))];
% %             KlinVoLin = II'\Ys(:,1);
% %             Klin = KlinVoLin(1:end-1);
% %             VoLin=KlinVoLin(end);
% %             clf;
% %             plot(Klin)
% %             title ('Composite kernel')
% %             saveas(1, sprintf('Composite kernel_M%d.pdf',EM.M))
            
            paramTrue = [];            
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% fit for each repetition
    %T = size(Iapps,1);
    for rep = 1:EXP.numRep
        disp(sprintf('Fit %dth course',rep))

        % read one repetition
        Iapp = Iapps(1:EM.T,rep);
        Y = Ys(1:EM.T,rep);
        if exist('Y2s','var')
            Y2 = Y2s(1:EM.T,rep);
        end
        if (EXP.TRUE)
            Iinj = Iinjs(1:EM.T,rep);
            V = Vms(1:EM.T,rep);
        end
        
        
        %% estimate composite kernel 
        disp('estimate composite kernel')
        tic
        II=stackCols(Iapp,EM.M*4,0);
        %II=stackCols(Iapps(:,rep),EM.M*4,0);  % should use longer data to calc kernel accurately for analysis
        II = [II; ones(1,size(II,2))];
        
        KlinVoLin = II'\Y;
        %KlinVoLin = II'\Ys(:,rep);
        Klin(:,rep,i) = KlinVoLin(1:end-1);
        VoLin(rep,i) = KlinVoLin(end);
        toc
        
        % estimate kernel from reference electrode!
        if exist('Y2s','var')
            KrefVoLin = II'\Y2;     % use all the samples!
            %KrefVoLin = II'\Y2s(:,rep);     % use all the samples!
            Kref(:,rep,i) = KrefVoLin(1:end-1);
            VoRef(rep,i) = KrefVoLin(end);
            
        end
%         if rep==1
%             clf;
%             plot(Klin(:,rep,i))
%             title ('linear kernel ')
%             
%             if exist('Kref','var')
%                 hold on
%                 plot(Kref(:,rep,i),'g--');
%                 legend ('composite','reference')
%             end
%             
%             
%             saveas(1, sprintf('%s_rep%d_composite_kernel_M%d.pdf',EXP.baseFilename,rep,EM.M))
%         end
            
        
        %%
        switch (EXP.METHOD)
            case 'AEC'
                disp('run AEC method')
                % run AEC from here !
                AEC.M = 50;
                if exist('Ke','var')
                    KeTrue = [Ke zeros(1,AEC.M-length(Ke))];
                end
                AEC.Mtot = AEC.M*4;
                AEC.Mtail = AEC.M+1;
                
                % 1. estimat composite kernel
                disp('estimate composite kernel')
                II=[stackCols(Iapp(:),AEC.Mtot,0); ones(1,length(Iapp))];
                K = II'\Y(:);

                % 2. estimate Ke from K
                try
                    [AEC.Ke AEC.Km] = electrodeKernel(K(1:end-1),AEC.Mtail);
                    AEC.Ke = AEC.Ke(1:AEC.M);
                catch exception
                    % statements
                    disp('error occured in AEC fitting')
                end
                

                %% filtering 
                U = filter(AEC.Ke, 1, Iapp);
                Vest = Y - U;

                if (EXP.TRUE)
                    mseKe(rep,i) = (AEC.Ke-KeTrue(:))'*(AEC.Ke-KeTrue(:))/AEC.M;    
                    
                %if(~isempty(V))
                    % not include initial Next points into consideration
                    % give time for steady state
                    %mseEMKeLin(rep,i) = (V(EM.Next+1:end)'-Vest(EM.Next+1:end)')*(V(EM.Next+1:end)-Vest(EM.Next+1:end))/length(V(EM.Next+1:end));
                    mse(rep,i) = (V(EM.Next+1:end)'-Vest(EM.Next+1:end)')*(V(EM.Next+1:end)-Vest(EM.Next+1:end))/length(V(EM.Next+1:end));
                end

                
                %% save AEC specific results 
                KeEsts(:,rep,i) = AEC.Ke;
                
            case {'sub','spike-single'}     % do single-tap fitting by EM     
                if ESTIMATE_SINGLE_TAP_EM
                    disp('init Ke from AEC')
                    
                    
                    %% RUN-AEC and use estimate Ke as a starting poit (July 1, 2012)
                    
                    if EM.INIT_KE_BY_AEC
                        try
                            [AEC.Ke AEC.Km] = electrodeKernel(Klin(:,rep,i),EM.M+1);
                            AEC.Ke = AEC.Ke(1:EM.M);

                            initKe = AEC.Ke;
                        catch exception
                            % statements
                            error('error occured in AEC fitting')
                            %initKe = [];
                        end

                        % filtering with AEC results 
                        U = filter(AEC.Ke, 1, Iapp);
                        VestByAEC = Y - U;
                        initVo = mean(VestByAEC);

                        % store AEC results 
                        VestsByAEC(:,rep,i) = VestByAEC;

                        % calc MSE for AEC
                        if (EXP.TRUE)
                            mseByAEC(rep,i) = calcMse(VestByAEC(EM.Next+1:end),V(EM.Next+1:end),0);
                        end

                        if exist('Y2','var')
                            mseByAEC(rep,i) = calcMse(VestByAEC(EM.Next+1:end),Y2(EM.Next+1:end),1);
                        end
                    else
                        initKe = [];
                        initVo = mean(Y);
                    end
                    
                    
                    
                    disp('run single-tap EM')
                    
                    %% data specific initialization & call EM 
                    switch (EXP.DATA)

% %                         case {'LDS','LDS_S'}
% %                             initGuess.alpha = paramTrue.alpha;
% %                             initGuess.beta = paramTrue.beta;
% %                             initGuess.gamma = paramTrue.gamma;
% %                             initGuess.R = paramTrue.R;
% %                             initGuess.Q = paramTrue.Q;
% %                             initGuess.Xo = -60;
% %                             initGuess.Po = 5;
% % 
% %                             paramEst= fit_EM(Iapp(:),Y(:)', dt, initGuess, EM, param);

                        case 'HH'
            %                 
            %                 initGuess.alpha = .9;
            %                 initGuess.beta = .1/5.55;
            %                 initGuess.gamma = 0.3*0.1*10.6;
            %                 initGuess.R = R;
            %                 initGuess.Q = 5;
            %                 initGuess.Xo = -60;
            %                 initGuess.Po = 5;

                            %% call the new EM routine
                            paramLow =[0.9, 1e-3, initVo-5, 2, R*2];
                            paramHigh =[0.99, 1e-1, initVo+5, 10, R*10];

                            initPrior = 10;


                        case 'CellB'
                            paramLow =[0.9, 1e-3, mean(Y)-5, 1e-2, 1e-2];
                            paramHigh =[0.99, 1e-1, mean(Y)+5, 1e-1, 1e-1];

                            initPrior = 1;


                        case 'Cell2'
            %                 initGuess.alpha = .99;
            %                 initGuess.beta = dt/sum(Klin);
            %                 initGuess.gamma = -0.1*dt*67;
                            initGuess.alpha = 0.99;
                            initGuess.beta = 0.01;
                            initGuess.gamma = -0.15;
                            initGuess.R = 10;
                            initGuess.Q = 10;
                            initGuess.Xo = mean(Y);
                            initGuess.Po = 10;

                            [paramEst Klin] = fit_EM(Iapp(:),Y(:)', initGuess, EM);                 
                        otherwise
                            error ('init is not done approprieately')
                                
                            
                    end

                    % store init and estimate for individual trials (Apr 30 2012)
                    [paramEst KeEst LL LLs initParams initKes allParams allKes] = estimParamEM(Y, Iapp, paramLow, paramHigh, initPrior, initKe, EM.M, EM.max_iter, EM.num_trial);
                    LLtrials = max(LLs);
                    itrTrials = sum(~isnan(LLs));



                    %% store EM-specific estimation results
                    maxItrs(rep,i) = max(itrTrials);
                    paramEsts(:,rep,i) = paramEst;
                    KeEsts(:,rep,i) = KeEst;
                    LLfinals(rep,i) = LL;
                    [~,maxTrial] = max(max(LLs));



                    % filtered input
                    U = filter(KeEst, 1, Iapp);
                    %% subtract filted Ve
                    VestEMKeLin = Y - U;

                    if (EXP.TRUE)
                    %if(~isempty(V))
                        % not include initial Next points into consideration
                        % give time for steady state
                        %mseEMKeLin(rep,i) = (V(EM.Next+1:end)'-Vest(EM.Next+1:end)')*(V(EM.Next+1:end)-Vest(EM.Next+1:end))/length(V(EM.Next+1:end));
                        mseEMKeLin(rep,i) = calcMse(VestEMKeLin(EM.Next+1:end),V(EM.Next+1:end),0);
                    end
                    
                    if exist('Y2','var')
                        mseEMKeLin(rep,i) = calcMse(VestEMKeLin(EM.Next+1:end),Y2(EM.Next+1:end),1);
                    end


                    %% Smoothe with estimated parameters
                    %[Xs, Ps, Pcs Xf Pf] = kalman_smth(Y'-paramEst(3), U', paramEst(1), paramEst(2), 1, 1, paramEst(4), paramEst(5), 0, 10);
                    [Xs, Ps] = kalman_smth(Y'-paramEst(3), U', paramEst(1), paramEst(2), 1, 1, paramEst(4), paramEst(5), 0, initPrior);
                    Vest = Xs' + paramEst(3);

                    % calc MSE
                    if (EXP.TRUE)
                        % MSE of Ke
                        mseKe(rep,i) = (KeEst-Ke)*(KeEst-Ke)'/length(Ke);
                        
                        % MSE of Vm
                        % not include initial Next points into consideration
                        % give time for steady state
                        mse(rep,i) = calcMse(Vest(EM.Next+1:end),V(EM.Next+1:end),0);
                    end
                    
                    if exist('Y2','var')
                        mse(rep,i) = calcMse(Vest(EM.Next+1:end),Y2(EM.Next+1:end),1);
                    end

                end % end of if ESTIMATE_SINGLE_TAP_EM

        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            case 'spike'    % for multi-tap fitting 
                %% 
                disp('read estimated params from single-tap fitting')
                paramEst = singleTap.paramEsts(:,rep,i);
                KeEst = singleTap.KeEsts(:,rep,i);
%                 LL = singleTap.LLfinals(rep,i)
                initPrior = 10;     % ??

                % filtered input
                U = filter(KeEst, 1, Iapp);

                %% Smoothe with single-tap results 
                % paramEst = [apha beta Vo Q R]
                disp('smoothing with single-tap results')
                [Xs, Ps] = kalman_smth(Y'-paramEst(3), U', paramEst(1), paramEst(2), 1, 1, paramEst(4), paramEst(5), 0, initPrior);
                Vest = Xs' + paramEst(3);
        
                % calc MSE
                if (EXP.TRUE)
                %if(~isempty(V))
                    % not include initial Next points into consideration
                    % give time for steady state
                    mseSingleTap(rep,i) = calcMse(Vest(EM.Next+1:end),V(EM.Next+1:end),0);
                    
                    %% estimate "ground truth" params by regression from true V
                    disp('ground truth from true V')
                    betaHat = paramEst(2);
                    MM = fliplr(makeStimRows([0;V(1:end-1)-betaHat*U(1:end-1)],dimX)); % make design matrix for regressing Y on its history
                    %aT2 = M\Y  % solve regression as initial values for alphas
                    %Q2 = var(Y-Vo2-MM*aT2)  % use residual variance to initialize estimate for Q
                    MM = [MM ones(size(MM,1),1)];  % mean is also estimated!
                    aT3Vo=MM\(V-betaHat*U);  
                    aT3 = aT3Vo(1:end-1)
                    disp('corresponding roots are')
                    disp(roots([1;-aT3]))
                    Vo3 = aT3Vo(end)
                    Q3 = min(var(V-betaHat*U-MM*aT3Vo),10)  % use residual variance to initialize estimate for Q 
                    Q3 = max(Q3,1e-9);
                    
                    % store 
                    linRegMultiTrue(:,rep) = [aT3; Vo3; Q3];
                    
                    paramEstTrueV(:,rep,i) = [aT3; betaHat; Vo3; Q3; paramEst(end)];
                    %% Smoothe with true param results 
                    disp('smoothing with true param')
                    
                    VestTrueParam = smoothMulti(aT3, [betaHat Vo3 Q3 paramEst(end)], initPrior, Y, U); %% replace with  a function
                    mseTrueV(rep,i) =  calcMse(VestTrueParam(EM.Next+1:end),V(EM.Next+1:end),0);
                end
                
                if exist('Y2','var')
                    mseSingleTap(rep,i) = calcMse(Vest(EM.Next+1:end),Y2(EM.Next+1:end),1);
                end
                
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
                Q2 = residue-paramEst(end);  % use residual variance to initialize estimate for Q 
                if Q2>10
                    fprintf('[WARNING] Q is too large');
                    Q2 = min(Q2,10);
                end
                Q2 = max(Q2,1e-9);

                % store 
                linRegMulti(:,rep) = [aT2; Vo2; Q2];
                paramEstLinRegMulti(:,rep,i) = linRegMulti(:,rep);
                mseLinRegMulti(:,rep,i) = residue;
                    
                
                %% range of initial params
                AtLow = aT2-abs(aT2*0.1);
                AtHigh = aT2+abs(aT2*0.1);

                % use beta and R from linear regime
                initVo = mean(Y-U)+Vo2;
                ParamLow =[paramEst(2), initVo+5, Q2/2, paramEst(end)];
                ParamHigh =[paramEst(2), initVo-5, Q2*2, paramEst(end)];  
                %initPrior = 10;
                
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Run multi-tap EM
                if ESTIMATE_MULTI_TAP_EM
                    %% Re-fit with multi-tap EM (KeEst is fixed)
                    disp('Re-fit with multi-tap EM (beta, KeEst, R are fixed)')

                    %% Now, call EM to estimate parameters
                    % LLs are implemented Apr. 25, 2012
                    % inital and final estiamtes for all trials are added
                    % Apr. 29, 2012
                    [AtEst, paramEstMulti, LL, LLs, initAt, initParams, allAt, allParams] = estimParamEM_multi(Y, U, AtLow, AtHigh, ParamLow, ParamHigh, initPrior, EM.max_iter, EM.num_trial);
                    LLtrials = max(LLs);
                    itrTrials = sum(~isnan(LLs))
                    
%                     % ultimately I'll store all LLs
%                     % but for now let's ignore
%                     [AtEst, paramEstMulti, LL] = estimParamEM_multi(Y, U, AtLow, AtHigh, ParamLow, ParamHigh, initPrior, EM.max_iter, EM.num_trial)
%                     LLs = [];
                    
                    
                    % store results
                    AtEsts(:,rep,i) = AtEst;
                    paramEstMultis(:,rep,i) = paramEstMulti;
                    LLfinals(rep,i) = LL;
                
                    %% Smoothe with estimated parameters
                    %disp('smoothing with multi-tap results')
                    Vest = smoothMulti(AtEst, paramEstMulti, initPrior, Y, U); %% replace with  a function
                    
                    % calc MSE
                    if (EXP.TRUE)
                    %if(~isempty(V))
                        % not include initial Next points into consideration
                        % give time for steady state
                        mseMultiTap(rep,i) = calcMse(Vest(EM.Next+1:end),V(EM.Next+1:end),0);
                    end
                    
                    if exist('Y2','var')
                        mseMultiTap(rep,i) = calcMse(Vest(EM.Next+1:end),Y2(EM.Next+1:end),1);
                    end
                    
                    
                                      
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Run multi-tap GD
                if ESTIMATE_MULTI_TAP_GD
                    
                    %% prepare initial
                    if ESTIMATE_MULTI_TAP_EM   
                        % use exactly the same init as EM
                    else
                        %draw samples from the given range of init
                        initAt = repmat(AtLow,1,GD.num_trial) + repmat(AtHigh-AtLow,1,GD.num_trial).*rand(dimX,GD.num_trial);
                        initParams = repmat(ParamLow(:),1,GD.num_trial) + repmat(ParamHigh(:)-ParamLow(:),1,GD.num_trial).*rand(length(ParamLow),GD.num_trial);
                        
                    end
                        
                    %% estimtae multiple times
                    prshat = nan*ones(dimX+2,GD.num_trial);
                    GD_LLs = nan*ones(GD.num_trial,1);
                    for trial = 1:GD.num_trial
                        
                        % Alter initial values so roots are all abs<1
                        aT0 = initAt(:,trial);
                        rts = roots([1;-aT0]);
                        if any(abs(rts)>=1)
                            rts = rts./((max(abs(rts))+.01));
                        end
                        aT1 = poly(rts);
                        aT1 = -aT1(2:end)';



                        % embed these params in the parameter vector
                        %prs0 =  [aT1; log(1./Q2); log(1./initParams(4,trial));initParams(2,trial);log(initParams(1,trial))];
                        prs0 =  [aT1; -log(initParams(3,trial)); initParams(2,trial)]; % aT -log(Q2) Vo
                        % changed to fix beta, R
                        % input arguments are now        prs,y,ntaps,II,R,beta)

                        betaFixed = initParams(1,trial);
                        Rfixed = initParams(4,trial);
                        lFunc = @(prs)(neglogmargLi_LDSvec_DynamPrsOnly(prs,Y(:),dimX,U(:),Rfixed,betaFixed));   % updated to fix beta and r
                        
                        opts = optimset('display', 'iter','largescale','off','maxiter',2000,'maxfuneval',25000,'tolfun',1e-16,'tolx',1e-16);
                        [prs negLL] = fminunc(lFunc,prs0, opts); % compute estimate: aT -log(Q) muy
                

                        % store estimates and LL
                        prshat(:,trial) = prs;
                        GD_LLs(trial) = -negLL;
                        
%                         idxNan = GD.LLs < -1e24;
%                         GD.LLs(idxNan) = NaN;
            
                    end
                    
                    %% choose estimate with max LL
                    [LL maxTrialGD] = max(GD_LLs);
                    
                    
                    % store results
                    GD.AtEst(:,rep,i) = prshat(1:dimX,maxTrialGD);

                    qHat = exp(-prshat(dimX+1,maxTrialGD));
                    mu_yHat = prshat(dimX+2,maxTrialGD);
                    GD.paramEst(:,rep,i) = [betaFixed; mu_yHat; qHat; Rfixed];
                    GD.LLfinals(rep,i) = LL;
                    
                    
                    %% Smoothe with estimated parameters
                    %disp('smoothing with multi-tap results')
                    VestGD = smoothMulti(prshat(1:dimX,maxTrialGD), [betaFixed; mu_yHat; qHat; Rfixed], initPrior, Y, U); %% replace with  a function
                    
                    % calc MSE
                    if (EXP.TRUE)
                    %if(~isempty(V))
                        % not include initial Next points into consideration
                        % give time for steady state
                        %mseMultiTap(rep,i) = (V(EM.Next+1:end)'-VestGD(EM.Next+1:end)')*(V(EM.Next+1:end)-VestGD(EM.Next+1:end))/length(V(EM.Next+1:end));
                        GD.mseMultiTap(rep,i) = calcMse(VestGD(EM.Next+1:end),V(EM.Next+1:end),0);
                    end
                    
                    if exist('Y2','var')
                        GD.mseMultiTap(rep,i) = calcMse(VestGD(EM.Next+1:end),Y2(EM.Next+1:end),1);
                    end
                    
                    
                end
            case 'post-lin'                %% POST PROCESSING: EXP WITH LINEAR KERNEL (MAR. 22, 2012)

                %% read estimated params
                paramEst = singleTap.paramEsts(:,rep,i)
                KeEst = singleTap.KeEsts(:,rep,i)
                LL = singleTap.LLfinals(rep,i)

                % filtered input
                U = filter(KeEst, 1, Iapp);

                %% Smoothe with single-tap results 
                [Xs, Ps] = kalman_smth(Y'-paramEst(3), U', paramEst(1), paramEst(2), 1, 1, paramEst(4), paramEst(5), 0, 10);
                Vest = Xs' + paramEst(3);
                
                if (EXP.TRUE)
                    %if(~isempty(V))
                        % not include initial Next points into consideration
                        % give time for steady state
                        mseEM(rep,i) = (V(EM.Next+1:end)'-Vest(EM.Next+1:end)')*(V(EM.Next+1:end)-Vest(EM.Next+1:end))/length(V(EM.Next+1:end));
                end
                    
                %% subtract filted Ve
                VestEMKeLin = Y - U;
        
                if (EXP.TRUE)
                %if(~isempty(V))
                    % not include initial Next points into consideration
                    % give time for steady state
                    %mseEMKeLin(rep,i) = (V(EM.Next+1:end)'-Vest(EM.Next+1:end)')*(V(EM.Next+1:end)-Vest(EM.Next+1:end))/length(V(EM.Next+1:end));
                    mseEMKeLin(rep,i) = (V(EM.Next+1:end)'-VestEMKeLin(EM.Next+1:end)')*(V(EM.Next+1:end)-VestEMKeLin(EM.Next+1:end))/length(V(EM.Next+1:end));
                end

                % check difference between them.
                % (Vest-VestEMKeLin)'*(Vest-VestEMKeLin)/length(Vest)

% %                 %% for debug
% %                 clf;
% %                 subplot(211)
% %                 plot([KeEst(:) Ke(:)])
% %                 subplot(212)
% %                 plot([VestEMKeLin V])
% %                 title (sprintf('mse = %.3f', mseEMKeLin(rep,i)))
                
%                 %% just linear kernel 
%                 II=stackCols(Iapp,EM.M,0);
%                 II = [II; ones(1,length(Y))];
%                 KlinVo = II'\Y;
%                 %Klin = Klin-min(Klin);
%                 % remove DC bias
%                 KlinVo = KlinVo - [min(KlinVo(1:end-1))*ones(EM.M,1); 0]
%                 VestLinKe = transpose(KlinVo'*II);
%                 
%                 
%                 
%                 if (EXP.TRUE)
%                 %if(~isempty(V))
%                     % not include initial Next points into consideration
%                     % give time for steady state
%                     %mseEMKeLin(rep,i) = (V(EM.Next+1:end)'-Vest(EM.Next+1:end)')*(V(EM.Next+1:end)-Vest(EM.Next+1:end))/length(V(EM.Next+1:end));
%                     mseLinKe(rep,i) = (V(EM.Next+1:end)'-VestLinKe(EM.Next+1:end)')*(V(EM.Next+1:end)-VestLinKe(EM.Next+1:end))/length(V(EM.Next+1:end));
%                 end
%                 
%                 %%
%                 clf;
%                 subplot(211)
%                 plot([KlinVo(1:end-1) Ke(:)])
%                 subplot(212)
%                 plot([VestLinKe V])
                
                
                % check difference between them.
                % (Vest-VestEMKeLin)'*(Vest-VestEMKeLin)/length(Vest)
                
        end
        
        
        

        
        %% plot and save fitting results for the first rep
        if rep == 1
        %if 1==1    
            %% 
            %idx = 1:min(2000,EM.T);
            %idx = 1:EM.T;
            idx = (EM.Next+1):EM.T;
            clf;
            numRow=7;
            subplot(numRow,1,1);
            plot(dt*(idx-EM.Next),Iapp(idx),'r');hold on
            if (EXP.TRUE)
                plot(dt*(idx-EM.Next),Iinj(idx));
                legend('I_{app}','I_{inj}'); legend boxoff
            end
            xlabel('ms')
            ylabel('uA/cm^2')
            title ('applied current')
            box off

            subplot (numRow,1,2)
            plot (dt*(idx-EM.Next),Y(idx))
            if exist('Y2s','var')    
                hold on;
                plot (dt*(idx-EM.Next),Y2(idx),'g')
                legend ('measure', 'reference');
            end
            title ('measured voltage (V_{rec})')
            box off
            ylabel('mV')
    %         if (EXP.DATA(1:4) == 'Cell')
    %             set(gca,'ylim', [-80 -60])
    %         end


            subplot(numRow,1,3);
            if exist('Vest','var')
                plot(dt*(idx-EM.Next),Vest(idx));  hold on
            end
            %plot(dt*(idx-EM.Next),Y(idx)-U(idx),'g');
            
            if exist('V','var')
                plot(dt*(idx-EM.Next),V(idx),'r');
                legend('estimated','true'); legend boxoff
            elseif exist('Y2s','var')
                plot(dt*(idx-EM.Next),Y2(idx),'g');
                legend('estimated','reference'); legend boxoff
            end
                
            % plot title with MSE
            switch EXP.METHOD
                case {'spike'}
%                     if ~ESTIMATE_MULTI_TAP_EM
%                         title(sprintf('Membrane voltage (V_m) MSE=%.2f',mseSingleTap(rep,i)));
%                     else
                        title(sprintf('Membrane voltage (V_m) MSE=%.2f (single) MSE=%.2f (multi)',mseSingleTap(rep,i),mseMultiTap(rep,i)));
%                     end
                otherwise
                    title(sprintf('Membrane voltage (V_m) MSE=%.2f',mse(rep,i)));
            end
                        
                    
            
            
            xlabel('ms')
            ylabel('mV')
    %         if (EXP.DATA(1:4) == 'Cell')
    %             set(gca,'ylim', [-80 -60])
    %         end
            box off

            
            %if (ESTIMATE_MULTI_TAP_EM & EXP.MODE~=0)      % 
            if exist('LLs','var') 
                if ~isempty(LLs)
                    subplot(numRow,2,7)
                    plot(LLs)
                    set(gca,'xlim', [1 max(itrTrials)])
                    title ('LL')
                    xlabel ('iteration')
                    box off

                    subplot(numRow,2,8)
                    plot(2:(size(LLs,1)),diff(LLs))
                    set(gca,'xlim', [1 max(itrTrials)])
                    title ('diff(LL)')
                    xlabel ('iteration')
                    box off


                    subplot(numRow,2,9)
                    %hist(LLtrials)
                    plot (LLtrials, 'o--')
                    
                    % high light max LL
                    hold on;
                    [LLmax maxTrial] = max(LLtrials);
                    plot(maxTrial,LLmax,'ro')
                    
                    title ('final LLs')
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
                end
            end


            % plot init & estimates params
            switch EXP.METHOD
                case {'sub','spike-single'}
                    % initial and estimate for individual trials (Apr 30
                    % 2012)
                    subplot(numRow,6,31);
                    plot(initKes','g+:'); hold on
                    plot(allKes','bo:'); hold on
                    xlabel('trial')
                    title ('Ke')
                    box off

                    subplot(numRow,6,32);
                    plot(initParams(1,:)','g+:'); hold on
                    plot(allParams(1,:)','bo:');
                    plot(maxTrial,allParams(1,maxTrial), 'ro')
                    xlabel('trial')
                    title ('\alpha')
                    box off

                    subplot(numRow,6,33);
                    plot(initParams(2,:)','g+:'); hold on
                    plot(allParams(2,:)','bo:'); hold on
                    plot(maxTrial,allParams(2,maxTrial), 'ro')
                    xlabel('trial')
                    title ('\beta')
                    box off
                    set(gca,'yscale','log')

                    subplot(numRow,6,34);
                    plot(initParams(3,:)','g+:'); hold on
                    plot(allParams(3,:)','bo:'); hold on
                    plot(maxTrial,allParams(3,maxTrial), 'ro')
                    xlabel('trial')
                    title ('V_{rev}')
                    box off
                    
                    subplot(numRow,6,35);
                    plot(initParams(4,:)','g+:'); hold on
                    plot(allParams(4,:)','bo:'); hold on
                    plot(maxTrial,allParams(4,maxTrial), 'ro')
                    xlabel('trial')
                    title ('Q')
                    set(gca,'yscale','log')
                    box off

                    subplot(numRow,6,36);
                    plot(initParams(5,:)','g+:'); hold on
                    plot(allParams(5,:)','bo:'); hold on
                    plot(maxTrial,allParams(5,maxTrial), 'ro')
                    xlabel('trial')
                    title ('R')
                    set(gca,'yscale','log')
                    box off
                    legend ('initial','estimate'); legend boxoff
                    
                        
                        
                    
                    % final estimate
                    
                    subplot(numRow,2,13);       % KeEst
                    plot (Klin(:,rep,i), 'k--'); hold on
                    plot (initKe,'g-.')
                    if exist('Ke','var')
                        plot([KeEst(:) Ke(:) ]);
                        title(sprintf('Estimated kernel (MSE=%.2e)',mseKe(rep,i)));
                    else
                        plot([ KeEst(:)]);
                        title(sprintf('Estimated kernel'));
                    end
                    set(gca,'xlim',[1 EM.M])
                    ylim = get(gca,'ylim');
                    set(gca,'ylim',[0 ylim(end)])
                    
                    legend ('composite', 'init', 'estim','true'); legend boxoff
                    box off
                    
                    subplot(numRow,2,14);       % other params
                    plot(paramEst);
                    title ('estimated params')
                    box off

                    
                    box off

                case 'spike'
                    if (ESTIMATE_MULTI_TAP_EM)
                        subplot(numRow,5,26);
                        plot(initAt','g+:'); hold on
                        plot(allAt','bo:'); hold on
                        if (EXP.TRUE)
                            plot(repmat(linRegMultiTrue(1:dimX,rep),1,EM.num_trial)','r--');    %regression with true V
                        end
                        
                        xlabel('trial')
                        title ('aT')
                        box off
                        
                        subplot(numRow,5,27);
                        plot(initParams(1,:)','g+:'); hold on
                        plot(allParams(1,:)','bo:'); hold on
                        xlabel('trial')
                        title ('\beta')
                        box off
                        set(gca,'yscale','log')
                        
                        subplot(numRow,5,28);
                        plot(initParams(2,:)','g+:'); hold on
                        plot(allParams(2,:)','bo:'); hold on
                        if (EXP.TRUE)
                            plot(repmat(linRegMultiTrue(dimX+1,rep),1,EM.num_trial)','r--'); %regression with true V
                        end
                        xlabel('trial')
                        title ('V_{rev}')
                        box off
                        
                        subplot(numRow,5,29);
                        plot(initParams(3,:)','g+:'); hold on
                        plot(allParams(3,:)','bo:'); hold on
                        if (EXP.TRUE)
                            plot(repmat(linRegMultiTrue(dimX+2,rep),1,EM.num_trial)','r--'); %regression with true V
                        end
                        xlabel('trial')
                        title ('Q')
                        set(gca,'yscale','log')
                        box off
                        
                        subplot(numRow,5,30);
                        plot(initParams(4,:)','g+:'); hold on
                        plot(allParams(4,:)','bo:'); hold on
                        xlabel('trial')
                        title ('R')
                        set(gca,'yscale','log')
                        box off
                        legend ('initial','estimate'); legend boxoff
                        
                        
                        subplot(numRow,2,13);
                        plot(AtEst); hold on
                        plot([paramEst(1) zeros(1,dimX-1)],'kx--')
                        title ('estimated aT with the largest LL')
                        box off
                        legend ('multi','single'); legend boxoff
                        
                        subplot(numRow,2,14);
                        plot(paramEstMulti); hold on
                        plot(paramEst(2:end),'kx--')
                        title ('estimated param with the largest LL')
                        box off
                    end
                    
                    
                    %% plot multi-tap estimates by GD
                    if (ESTIMATE_MULTI_TAP_GD)
                        LINE_TYPE_GD = 'ks:'
                        subplot(numRow,5,26);   % aT for all trial
                        plot(initAt','g+:'); hold on
                        plot(prshat(1:dimX,:)',LINE_TYPE_GD); hold on
%                         if (EXP.TRUE)
%                             plot(repmat(linRegMultiTrue(1:dimX,rep),1,GD.num_trial)','r--');    %regression with true V
%                         end
%                         
                        xlabel('trial')
                        title ('aT')
                        box off
                        
                        subplot(numRow,5,27); % beta for all trial (no need)
                        plot(initParams(1,:)','g+:'); hold on
                        %plot(allParams(1,:)','bo:'); hold on
                        xlabel('trial')
                        title ('\beta')
                        box off
                        set(gca,'yscale','log')
                        
                        subplot(numRow,5,28);
                        plot(initParams(2,:)','g+:'); hold on
                        plot(prshat(dimX+2,:)',LINE_TYPE_GD); hold on
                        if (EXP.TRUE)
                            plot(repmat(linRegMultiTrue(dimX+1,rep),1,GD.num_trial)','r--'); %regression with true V
                        end
                        xlabel('trial')
                        title ('V_{rev}')
                        box off
                        
                        subplot(numRow,5,29);
                        plot(initParams(3,:)','g+:'); hold on
                        plot(prshat(dimX+1,:)',LINE_TYPE_GD); hold on
                        if (EXP.TRUE)
                            plot(repmat(linRegMultiTrue(dimX+2,rep),1,GD.num_trial)','r--'); %regression with true V
                        end
                        xlabel('trial')
                        title ('Q')
                        set(gca,'yscale','log')
                        box off
                        
                        subplot(numRow,5,30);
                        plot(initParams(4,:)','g+:'); hold on
                        %plot(allParams(4,:)','bo:'); hold on
                        xlabel('trial')
                        title ('R')
                        set(gca,'yscale','log')
                        box off
                        legend ('initial','estimate'); legend boxoff
                        
                        
                        subplot(numRow,2,13);
                        plot(GD.AtEst(1:dimX,rep,i),LINE_TYPE_GD); hold on
                        %plot([paramEst(1) zeros(1,dimX-1)],'kx--')
                        title ('estimated aT with the largest LL')
                        box off
                        legend ('multi','single'); legend boxoff
                        
                        subplot(numRow,2,14);
                        plot(GD.paramEst(dimX+1:end,rep,i),LINE_TYPE_GD); hold on
                        title ('estimated param with the largest LL')
                        box off
                    end
            end
                
            set(gcf, 'paperposition', [0 0 12, 20]) 
            set(gcf, 'papersize', [12 20]) 
            
            switch EXP.METHOD
                case 'spike'
                    % save figure
                    filename = sprintf('%s%s_R%.2f_dim%d_N%d_itr%d_trial%d_rep%03d.pdf',EXP.DATA,EXP.METHOD,R,dimX,EM.T,EM.max_iter,EM.num_trial,rep)                    
                    saveas(1, filename)
                    
                    filename = sprintf('%s%s_R%.2f_dim%d_N%d_itr%d_trial%d_rep%03d.mat',EXP.DATA,EXP.METHOD,R,dimX,EM.T,EM.max_iter,EM.num_trial,rep)
                    
                otherwise 
                    % save figure
                    filename = sprintf('%s%s_R%.2f_N%d_itr%d_trial%d_rep%03d.pdf',EXP.DATA,EXP.METHOD,R,EM.T,EM.max_iter,EM.num_trial,rep)                    
                    saveas(1, filename)
                    
                    filename = sprintf('%s%s_R%.2f_N%d_itr%d_trial%d_rep%03d.mat',EXP.DATA,EXP.METHOD,R,EM.T,EM.max_iter,EM.num_trial,rep)
                    
            end
            
            
%            if ESTIMATE_MULTI_TAP_EM | ESTIMATE_MULTI_TAP_GD
                
% %             % save only results for current rep 
% %             if exist('paramEstMulti','var')
% %                 if exist('LLs','var')
% %                     switch EXP.METHOD
% %                         case 'spike'
% %                             % to save all initial and intermediate results 
% %                             if exist('linRegMultiTrue','var')
% %                             
% %                                 save(filename,'paramEst','KeEst','paramEstMulti','LLs','dimX','initParams','initAt','allParams','allAt','linRegMulti','linRegMultiTrue')
% %                             else
% %                                 save(filename,'paramEst','KeEst','paramEstMulti','LLs','dimX','initParams','initAt','allParams','allAt','linRegMulti')
% %                             end
% %                             
% %                         otherwise
% %                             save(filename,'paramEst','KeEst','paramEstMulti','LLs','dimX','initParams','initKes','allParams','allKes')
% %                     end
% %                 else
% %                     save(filename,'paramEst','KeEst','paramEstMulti')
% %                 end
% %             else
% %                 save(filename,'paramEst','KeEst')
% %             end
            
        end  % end of plotting for rep==0
        
        
        %% store estimation results
        Vests(:,rep,i) = Vest;
        
        % store estimated params
        %paramEsts(rep,i) = paramEst;
        % into matrix
        %params(:,rep,i) = [paramEst.alpha, paramEst.beta, paramEst.gamma, paramEst.Ke, paramEst.Q, paramEst.R, paramEst.Xo, paramEst.Po]';
        
    end


    timeDuration(i) = toc
end


%% save results
if exist('VestsByAEC','var')
    save VestsByAEC
end

%% Save linear kernels into a mat file 

filename = sprintf('%s_linear_kernel_rep1-%d.mat',EXP.DATA, EXP.numRep)
if exist('Kref','var')
    save(filename,'EM','EXP','Rs','Klin','VoLin', 'Kref','VoRef')
else
    save(filename,'EM','EXP','Rs','Klin','VoLin')
end

%%
switch EXP.METHOD
    case {'sub','spike-single'}
        filename = sprintf('fit_%s%s_rep1-%d_N%d_itr%d_trial%d.mat',EXP.DATA,EXP.METHOD,EXP.numRep,EM.T,EM.max_iter,EM.num_trial)
        save (filename, 'EXP', 'Qs', 'Rs', 'paramEsts', 'KeEsts', 'paramTrue', 'Vests', 'VestsByAEC', 'mse', 'mseByAEC', 'mseEMKeLin', 'mseKe', 'LLfinals', 'EM', 'maxItrs','timeDuration')
    case 'spike'
        filename = sprintf('fit_%s%s_rep1-%d_dim%d_N%d_itr%d_trial%d.mat',EXP.DATA,EXP.METHOD,EXP.numRep,dimX,EM.T,EM.max_iter,EM.num_trial)
        if ESTIMATE_MULTI_TAP_EM
            save (filename, 'EXP', 'Qs', 'Rs', 'paramEstMultis', 'AtEsts', 'paramTrue', 'Vests', 'mseSingleTap', 'mseMultiTap', 'LLfinals', 'EM', 'maxItrs', 'dimX','timeDuration')
        else
            filename = sprintf('fit_%s%s_rep1-%d_dim%d_N%d_itr%d_trial%d.mat',EXP.DATA,EXP.METHOD,EXP.numRep,dimX,EM.T,EM.max_iter,EM.num_trial)
            save (filename, 'EXP', 'Qs', 'Rs', 'paramTrue', 'mseSingleTap', 'dimX','timeDuration')
        end
        
        if ESTIMATE_MULTI_TAP_GD
            filename = sprintf('fit_%s%s_GD_rep1-%d_dim%d_N%d_itr%d_trial%d.mat',EXP.DATA,EXP.METHOD,EXP.numRep,dimX,EM.T,EM.max_iter,EM.num_trial)
            save (filename, 'EXP', 'Qs', 'Rs', 'GD', 'paramEstTrueV', 'paramEstLinRegMulti', 'dimX', 'timeDuration','Vests', 'mseSingleTap','mseTrueV','mseLinRegMulti','mseMultiTap')
        end
          
        
        
    case 'AEC'
        filename = sprintf('fit_%s%s_rep1-%d_N%d.mat',EXP.DATA,EXP.METHOD,EXP.numRep,EM.T)
        %save (filename, 'Qs', 'Rs', 'paramEsts', 'KeEsts', 'paramTrue', 'mse', 'mseKe', 'EM')
        save (filename, 'EXP', 'AEC','EM','Qs', 'Rs',  'paramTrue','KeEsts', 'Vests', 'mse', 'mseKe', 'timeDuration')
        
end


if exist('paramEstTrueV','var')
    filename = sprintf('paramEstTrueV_%s%s_rep1-%d_dim%d_N%d_itr%d_trial%d.mat',EXP.DATA,EXP.METHOD,EXP.numRep,dimX,EM.T,EM.max_iter,EM.num_trial)
    save (filename, 'Qs', 'Rs',  'paramEstTrueV', 'mse', 'mse', 'EM','timeDuration')

end
    


% %% plot LLs  %=> Doesn't have much meaning
% clf;
% if (size(LLs,1)~=1)
%     errorbar(Rs(1:numR), mean(LLs), std(LLs))
% else
%     semilogx(Rs, LLs,'+-')
% end
% ylabel('LL')
% xlabel('R')
% filename = sprintf('LL_%s_Mstep%d%d_T%d.pdf',EXP.DATA,EM.MstepMethod,EM.MstepConstraint,EM.T)
% saveas(1, filename)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot mse as function of R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load fit_HHsub_rep1-10_N5100_itr2000_trial32
% load fit_HHsub_rep1-10_N10100_itr2000_trial32

%load fit_HHspike_rep1-10_N5100_itr5000_trial32
%load fit_HHspike_rep1-10_N10100_itr5000_trial32


%


% linear: single tap
%load fit_HHsub_rep1-10_N5100_itr2000_trial32.mat
%load fit_HHsub_rep1-50_N5100_itr5000_trial32.mat
% 
% % spiking: single tap
%load fit_HHspike-single_rep1-10_N5100_itr5000_trial32.mat
% 
% % spiking: multi-tap
% load fit_HHspike_rep1-10_dim1_N5100_itr5000_trial32.mat
%load fit_HHspike_rep1-10_dim2_N5100_itr5000_trial32.mat; dimX=2
% load fit_HHspike_rep1-10_dim3_N5100_itr5000_trial32.mat; dimX = 3
% load fit_HHspike_rep1-10_dim4_N5100_itr5000_trial32.mat; dimX=4


if (EXP.TRUE)

    clf;
    set(0,'DefaultAxesLineWidth', 1)
    
    
    LINE_TYPE = 'k'
    LINE_TYPE2 = '--x';
    LINE_TYPE_LIN_REG='xk:';
    LINE_TYPE_MULTI_EM = 'o-';
    LINE_TYPE_MULTI_GD = 'sg-';
    PLOT_REFERENCE = 1;
    SET_AXIS_LOG2 = 1;
    %LINE_WIDTH = 1.5;
    %% 
    LINE_WIDTH = 1;    
    BAR_WIDTH_INC =  1.2;
    
    numR=length(Rs)
    switch EXP.METHOD
        case {'sub','spike-single','AEC'}
            subplot(221)
            if (size(mseKe,1)~=1)
                if SET_AXIS_LOG2
                    %hE = errorbar(log2(Rs(1:numR)), mean(log2(mseKe)), std(log2(mseKe)),LINE_TYPE,'lineWidth',LINE_WIDTH); hold on
                    hE = errorbar(log2(Rs(1:numR)), mean(mseKe), std(mseKe),LINE_TYPE,'lineWidth',LINE_WIDTH); hold on
                    xlabel('log_2(R)')
                    set(gca,'xlim', [-1.5 3.5])
                else 

                    hE = errorbar(Rs(1:numR), mean(mseKe), std(mseKe),LINE_TYPE,'lineWidth',LINE_WIDTH); hold on

                    xlabel('R')
                    set(gca,'xtick',Rs)
                end

            else
                plot(Rs, mseKe,'+-','lineWidth',LINE_WIDTH);                    
            end
            title('mse of kernel')

            box off



    %         set(gca,'xtick', Rs)
            set(gca,'ticklength', 4*[0.0100    0.0250])



            subplot(222)
            
            
        otherwise       % plot for spiking regime (multi-tap)
            subplot(221)
            
            
            
    end
    
    

    
    
    %% Plot MSE of estimated V
    clear legendStr;
    legendCnt=1;

    hold on
    if exist('mseSingleTap','var')
        errorbar(log2(Rs), mean(log2(mseSingleTap)),std(log2(mseSingleTap)),'c+--'); 
        legendStr{legendCnt} = 'single tap'; legendCnt=legendCnt+1;
    end
    if ESTIMATE_MULTI_TAP_EM & exist('mseMultiTap','var')
        errorbar(log2(Rs), mean(log2(mseMultiTap),1),std(log2(mseMultiTap),[],1),LINE_TYPE_MULTI_EM,'lineWidth',LINE_WIDTH)
        legendStr{legendCnt} = 'multi tap (EM)'; legendCnt=legendCnt+1;
    end
    if ESTIMATE_MULTI_TAP_GD & exist('GD.mseMultiTap','var')
        errorbar(log2(Rs), mean(log2(GD.mseMultiTap),1),std(log2(GD.mseMultiTap),[],1),LINE_TYPE_MULTI_GD,'lineWidth',LINE_WIDTH)
        legendStr{legendCnt} = 'multi tap (GD)'; legendCnt=legendCnt+1;
    end

    if exist('mseTrueV','var')
        errorbar(log2(Rs), mean(log2(mseTrueV)),std(log2(mseTrueV)),'r--'); 
        legendStr{legendCnt} = 'true V'; legendCnt=legendCnt+1;

    end
    %legend (legendStr); legend boxoff


%         %% plot mean square error of Ke without smoothing
%         if exist('mseEMKeLin','var')
%             
%             if SET_AXIS_LOG2
%                 hE = errorbar(log2(Rs(1:numR)), mean(log2(mseEMKeLin)), std(log2(mseEMKeLin)),LINE_TYPE2,'lineWidth',LINE_WIDTH); hold on
%             else
%                 hE = errorbar(Rs(1:numR), mean(mseEMKeLin), std(mseEMKeLin),LINE_TYPE2,'lineWidth',LINE_WIDTH); hold on
%                 %axis ([0 8.5 0 8.5])
%             set(gca,'xtick',Rs)
%             set(gca,'ytick',Rs)
%             end
%             legend ('smoothing','w/o smoothing','location','NW'); legend boxoff
%         end
% 

    
    if (PLOT_REFERENCE)
        if SET_AXIS_LOG2
            plot(log2(Rs), log2(Rs), 'k:','lineWidth',LINE_WIDTH)
        else
            plot(Rs, Rs, 'k:','lineWidth',LINE_WIDTH)
        end
    end
       
    set(gca,'ticklength', 4*[0.0100    0.0250])
    if SET_AXIS_LOG2
        title('log_2(mse of V_m) by smoothing')
    else
        title('mse of V_m')
    end
    
    box off
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot estimate Q and Rs 

    
    
    switch EXP.METHOD
        case {'sub','spike-single'}
            QestEMs = reshape(paramEsts(end-1,:,:),size(paramEsts,2),[])
            RestEMs = reshape(paramEsts(end,:,:),size(paramEsts,2),[])

            %clf;
            subplot(223)    % estimated Qs
            hold on
    
%             subplot(224)    % estimated Rs
%             if size(RestEMs,1)>1
%                 errorbar(log2(Rs), mean(log2(RestEMs)),std(log2(RestEMs)))
%             else 
%                 plot(log2(Rs), log2(RestEMs))
%             end
%             hold on;
%             plot(log2(Rs), log2(Rs), '--k')
%             box off
%             xlabel('log_2(R)')
%             title('log_2(R_{est})')            
            
        case 'spike'
            %clf;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% to plot Vest
            subplot(223); hold on            
            clear legendStr;
            legendCnt=1;
            
            VestTrueV = reshape(paramEstTrueV(3,:,:),EXP.numRep,[])
            VestLinRegMulti = reshape(paramEstLinRegMulti(end-1,:,:),EXP.numRep,[])

            errorbar(log2(Rs), mean(VestTrueV,1),std(VestLinRegMulti,[],1),'r--');
            legendStr{legendCnt} = 'true V'; legendCnt=legendCnt+1;
            
            errorbar(log2(Rs), mean(VestLinRegMulti,1),std(VestLinRegMulti,[],1),LINE_TYPE_LIN_REG);
            legendStr{legendCnt} = 'lin reg'; legendCnt=legendCnt+1;
            
            if (ESTIMATE_MULTI_TAP_EM)  % Qest by EM
                VestEMs = reshape(paramEstMultis(end-1,:,:),size(paramEstMultis,2),[])
                errorbar(log2(Rs), mean(VestEMs,1) , std(VestEMs,[],1),LINE_TYPE_MULTI_EM);
                legendStr{legendCnt} = 'multi-tap EM'; legendCnt=legendCnt+1;
                
% %                 RestEMs = reshape(paramEstMultis(end,:,:),size(paramEstMultis,2),[])
% %                 
% %                 errorbar(log2(Rs), mean(log2(QestEMs),1),std(log2(QestEMs),[],1)); hold on
            end
            if (ESTIMATE_MULTI_TAP_GD)  % Qest by GD
                VestGDs=reshape(GD.paramEst(2,:,:),EXP.numRep,[]);
                
                errorbar(log2(Rs), mean(VestGDs,1) , std(VestGDs,[],1),LINE_TYPE_MULTI_GD);
                legendStr{legendCnt} = 'multi-tap GD'; legendCnt=legendCnt+1;
                
            end

            legend (legendStr); legend boxoff
            
            box off
            xlabel('log_2(R)')
            title('V_{est}')
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % plot estimated Qs
            subplot(224)    % 
            hold on
    
            QestTrueV = reshape(paramEstTrueV(end-1,:,:),EXP.numRep,[])
            QestLinRegMulti = reshape(paramEstLinRegMulti(end,:,:),EXP.numRep,[])

            errorbar(log2(Rs), mean(log2(QestTrueV),1),std(log2(QestLinRegMulti),[],1),'r--');
            errorbar(log2(Rs), mean(log2(QestLinRegMulti),1),std(log2(QestLinRegMulti),[],1),LINE_TYPE_LIN_REG);
            
            if (ESTIMATE_MULTI_TAP_EM)  % Qest by EM
                QestEMs = reshape(paramEstMultis(end-1,:,:),size(paramEstMultis,2),[])
                RestEMs = reshape(paramEstMultis(end,:,:),size(paramEstMultis,2),[])
                
                errorbar(log2(Rs), mean(log2(QestEMs),1),std(log2(QestEMs),[],1),LINE_TYPE_MULTI_EM); hold on
            end
            if (ESTIMATE_MULTI_TAP_GD)  % Qest by GD
                QestGDs=reshape(GD.paramEst(3,:,:),EXP.numRep,[]);
                
                errorbar(log2(Rs), mean(log2(QestGDs),1) , std(log2(QestGDs),[],1),LINE_TYPE_MULTI_GD);
            end

            

            box off
            xlabel('log_2(R)')
            title('log_2(Q_{est})')
            
            
            
            
            


    end

    
    

    

    
    set(gcf, 'paperposition', [0 0 8 6]) 
    set(gcf, 'papersize', [8 6])
    switch EXP.METHOD
        case 'spike'            
            filename = sprintf('MSE_QR_%s%s_rep1-%d_MSE_dim%d_N%d_itr%d_trial%d.pdf',EXP.DATA,EXP.METHOD,EXP.numRep,dimX,EM.T,EM.max_iter,EM.num_trial)
        otherwise 
            filename = sprintf('MSE_QR_%s%s_rep1-%d_MSE_N%d_itr%d_trial%d.pdf',EXP.DATA,EXP.METHOD,EXP.numRep,EM.T,EM.max_iter,EM.num_trial)
    end
    saveas(1, filename)
end



timeDuration



return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE: compare MSEs of AEC vs EM (July, 9, 2012)
% 1. Linear regime 
% re-run after bug fix and code cleaning-up
% single file contains AEC, EM, and GD
clear


% HH neuron
%load ../2012-07-09_HH_AEC-DEC/fit_HHsub_rep1-50_N5100_itr100_trial32.mat

% real data 
load ../2012-07-09_CellB_AEC_DEC/fit_CellBsub_rep1-10_N5100_itr100_trial32.mat


clf
subplot(121)
plot(mse); hold on
plot(mseByAEC, '--')
switch EXP.DATA
    case 'HH'
        legend (sprintf('R = %.2f',Rs(1)),...
                sprintf('R = %.2f',Rs(2)),...
                sprintf('R = %.2f',Rs(3)),...
                sprintf('R = %.2f',Rs(4)),...
                sprintf('R = %.2f',Rs(5)))
            
    case 'CellB'
        legend (sprintf('\\sigma_{in} = %.2f',Rs(1)),...
                sprintf('\\sigma_{in} = %.2f',Rs(2)),...
                sprintf('\\sigma_{in} = %.2f',Rs(3)),...
                sprintf('\\sigma_{in} = %.2f',Rs(4)))
end
legend boxoff
title ('MSE')
box off
xlabel('repeat')

% subplot(223); hold on
% plot(dcDiff)
% plot(AEC_FIT.dcDiff, '--')
% title ('DC differene')
% box off


%idxToPlot=[1:4 6:10];
idxToPlot=1:size(mseByAEC,1);
subplot(122); hold on
errorbar(Rs, mean(mseByAEC(idxToPlot,:)),std(mseByAEC(idxToPlot,:)),'k--')
errorbar(Rs, mean(mse(idxToPlot,:)),std(mse(idxToPlot,:)),'k')
box off
legend ('AEC', 'AEC-EM'); legend boxoff

% subplot(224); hold on
% errorbar(Rs, mean(AEC_FIT.dcDiff),std(AEC_FIT.dcDiff),'k--')
% errorbar(Rs, mean(dcDiff),std(dcDiff),'k')
% title ('DC differene as function of \sigma_{in} (excluding outlier)')

switch EXP.DATA
    case 'HH'
        title ('MSE as function of R')
        xlabel('R')
        
        plot (Rs,Rs,':k')
            
    case 'CellB'
        title ('MSE as function of \sigma_{in}')
        xlabel('\sigma_{in}')
end






set(gcf, 'paperposition', [0 0 8 3]) 
set(gcf, 'papersize', [8 3]) 
filename = sprintf('%s_%s_MSE_AEC_vs_EM_itr%d.pdf',EXP.DATA,EXP.METHOD,EM.max_iter)
saveas(1,filename)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. for spiking-regime with multi-tap EM
% plotted for different tapsize 
clear
for dim = 1:4
    %filename = sprintf('2012-05-17_HH_nke4_additional/fit_HHspike_rep1-50_dim%d_N5100_itr500_trial32.mat',dim);
    %filename = sprintf('fit_HHspike_rep1-50_dim%d_N5100_itr50_trial32.mat',dim)
    
    
    % load EM results 
    filename = sprintf('../2012-07-09_CellB_AEC_DEC/fit_CellBspike_rep1-10_dim%d_N5100_itr100_trial32.mat',dim)
    load(filename, 'AtEsts', 'dimX','mseMultiTap', 'mseSingleTap','paramEstMultis', 'EXP', 'Rs', 'EM')
    
    mseMean(dim,:)=mean(mseMultiTap)
    mseStd(dim,:)=mean(mseMultiTap)
    
    mseSingleMean(dim,:)=mean(mseSingleTap)
    mseSingleStd(dim,:)=mean(mseSingleTap)
    
end

clf
numRow=length(Rs)/2;
numCol=2;
for i=1:length(Rs)
    subplot(numRow,numCol,i)
    errorbar(mseSingleMean(:,i),mseSingleStd(:,i),'--'); hold on 
    errorbar(mseMean(:,i),mseStd(:,i)); hold on 
    
    
    if i==1
        legend ('single-tap','multi-tap')
    end
    
    switch EXP.DATA
        case 'HH'
            title (sprintf('MSE with R = %.2f',Rs(i)))
            
        case 'CellB'
            title (sprintf('MSE with \\sigma_{in} = %.2f',Rs(i)))
    end


    
    box off
    xlabel('tap')
    
end

set(gcf, 'paperposition', [0 0 4*numCol 3*numRow]) 
set(gcf, 'papersize', [4*numCol 3*numRow]) 
filename = sprintf('%s_%s_MSE_AEC_vs_EM_itr%d.pdf',EXP.DATA,EXP.METHOD,EM.max_iter)
saveas(1,filename)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compare MSE after removing DC bias (July, 8, 2012)
% old code where MSE needed to be calculated again 
clear

load ../2012-07-08_AEC_DEC_itr50/fit_CellBsub_rep1-10_N5100_itr50_trial32.mat
%load ../2012-07-08_AEC_DEC_itr100/fit_CellBsub_rep1-10_N5100_itr100_trial32.mat
%load ../2012-07-08_AEC_DEC_itr1000/fit_CellBsub_rep1-10_N5100_itr1000_trial32.mat

size(Vests)
len=size(Vests,1)



%AEC_FIT=load('../2012-06-29_EM_vs_AEC_linear/fit_CellB_AEC_rep1-10_N5100.mat')% Vests is not stored...
AEC_FIT=load('../2012-07-08_AEC/fit_CellBAEC_rep1-10_N5100.mat')% Vests is not stored...


% read Y2s
for i=1:length(Rs)
    filename = sprintf('../DataCell1/Cell_04_15_2010_BD_n%.2f.mat',Rs(i))
    load(filename,'Y2s')
    if iscell(Y2s)
        Y2s=cell2mat(Y2s);
    end
    Y2s=Y2s(1:len,:);
    
    
    [mse(:,i), dcDiff(:,i)] = calcMse(Vests(:,:,i),Y2s,1);
    
    [AEC_FIT.mse(:,i), AEC_FIT.dcDiff(:,i)] = calcMse(AEC_FIT.Vests(:,:,i),Y2s,1);
    
   
end

clf
subplot(221)
plot(mse); hold on
plot(AEC_FIT.mse, '--')
legend (sprintf('\\sigma_{in} = %.2f',Rs(1)),...
        sprintf('\\sigma_{in} = %.2f',Rs(2)),...
        sprintf('\\sigma_{in} = %.2f',Rs(3)),...
        sprintf('\\sigma_{in} = %.2f',Rs(4)))
legend boxoff
title ('MSE')
box off

subplot(223); hold on
plot(dcDiff)
plot(AEC_FIT.dcDiff, '--')
title ('DC differene')
box off


idxToPlot=[1:4 6:10];
subplot(222); hold on
errorbar(Rs, mean(AEC_FIT.mse(idxToPlot,:)),std(AEC_FIT.mse(idxToPlot,:)),'k--')
errorbar(Rs, mean(mse(idxToPlot,:)),std(mse(idxToPlot,:)),'k')
title ('MSE as function of \sigma_{in} (excluding outlier)')
box off
legend ('AEC', 'AEC-EM'); legend boxoff

subplot(224); hold on
errorbar(Rs, mean(AEC_FIT.dcDiff),std(AEC_FIT.dcDiff),'k--')
errorbar(Rs, mean(dcDiff),std(dcDiff),'k')
title ('DC differene as function of \sigma_{in} (excluding outlier)')
xlabel('\sigma_{in}')
box off


set(gcf, 'paperposition', [0 0 15, 12]) 
set(gcf, 'papersize', [15 12]) 
filename = sprintf('%s_MSE_AEC_vs_EM_itr%d.pdf',EXP.DATA,EM.max_iter)
saveas(1,filename)







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyze linear kernel  (June 27, 2012)
% This results will be use as ground truth for estimation with  real data

clear
%load ../2012-07-08_AEC_DEC/CellB_linear_kernel_rep1-10.mat
%load ../2012-07-08_AEC_DEC_itr50/CellB_linear_kernel_rep1-10.mat
load ../2012-07-08_AEC_DEC_itr100/CellB_linear_kernel_rep1-10.mat

KlinMean = reshape(mean(Klin,2),[],length(Rs));
KlinStd = reshape(std(Klin,[],2),[],length(Rs));
KrefMean = reshape(mean(Kref,2),[],length(Rs));
KrefStd = reshape(std(Kref,[],2),[],length(Rs));
KdiffMean = reshape(mean(Klin-Kref,2),[],length(Rs));
KdiffStd = reshape(std(Klin-Kref,[],2),[],length(Rs));

% save CellB_linear_kernel_rep1-10.mat

%%
subplot(221); 
errorbar(KlinMean,KlinStd)
legend (sprintf('\\sigma_{in} = %.2f',Rs(1)),...
        sprintf('\\sigma_{in} = %.2f',Rs(2)),...
        sprintf('\\sigma_{in} = %.2f',Rs(3)),...
        sprintf('\\sigma_{in} = %.2f',Rs(4)))
legend boxoff
title('composite kernel')
set(gca,'xlim',[1 size(KlinMean,1)])
set(gca,'ylim',[0 3.5])
%set(gca,'yscale','log')
box off
subplot(222)
errorbar(KrefMean,KrefStd)
title('reference kernel')
set(gca,'xlim',[1 size(KlinMean,1)])
box off
subplot(223)
plot(KdiffMean)
%errorbar(KdiffMean,KdiffStd)
title('mean of difference')
set(gca,'xlim',[1 50])
set(gca,'ylim',[0 3.5])
box off
subplot(224)
plot(KdiffStd)
title('std of difference')
set(gca,'xlim',[1 50])
box off

set(gcf, 'paperposition', [0 0 15, 12]) 
set(gcf, 'papersize', [15 12]) 
filename = sprintf('%s_linear_kernels_compare_100.pdf',EXP.DATA)
saveas(1,filename)


%% compare with kernels estimated by EM 
% clear
% load CellB_linear_kernel_rep1-10.mat

EM_FIT=load('../2012-06-10_CellB/fit_CellBsub_rep1-10_N5100_itr50_trial32.mat')
size(EM_FIT.KeEsts)
Rs = EM_FIT.Rs;
DATA=EM_FIT.EM.DATA

KemMean = reshape(mean(EM_FIT.KeEsts,2),[],length(Rs));
KemStd = reshape(std(EM_FIT.KeEsts,[],2),[],length(Rs));



AEC_FIT=load('../2012-06-29_EM_vs_AEC_linear/fit_CellB_AEC_rep1-10_N5100.mat')
KaecMean = reshape(mean(AEC_FIT.KeEstsByAEC,2),[],length(Rs));
KaecStd = reshape(std(AEC_FIT.KeEstsByAEC,[],2),[],length(Rs));


%AEC_EM_FIT=load('../2012-07-08_AEC_DEC_itr50/fit_CellBsub_rep1-10_N5100_itr50_trial32.mat')
AEC_EM_FIT=load('../2012-07-08_AEC_DEC_itr100/fit_CellBsub_rep1-10_N5100_itr100_trial32.mat')
KaecEmMean = reshape(mean(AEC_EM_FIT.KeEsts,2),[],length(Rs));
KaecEmStd = reshape(std(AEC_EM_FIT.KeEsts,[],2),[],length(Rs));


YLIM_FOR_STD=0.2;

clf

% plot AEC
subplot(421)
plot(KaecMean)
%errorbar(KdiffMean,KdiffStd)
title('mean of KeEst by AEC')
set(gca,'xlim',[1 50])
set(gca,'ylim',[0 3.5])
box off
subplot(422)
plot(KaecStd)
title('std')
set(gca,'xlim',[1 50])
set(gca,'ylim',[0 YLIM_FOR_STD])
box off

% plot EM
subplot(423)
plot(KemMean); 
%errorbar(KdiffMean,KdiffStd)
title('mean of KeEst by EM')
set(gca,'xlim',[1 50])
set(gca,'ylim',[0 3.5])
box off

subplot(424)
plot(KemStd)
title('std')
set(gca,'xlim',[1 50])
set(gca,'ylim',[0 YLIM_FOR_STD])
box off


% plot EM
subplot(425)
plot(KaecEmMean); 
%errorbar(KdiffMean,KdiffStd)
title('mean of KeEst by EM initialized by AEC')
set(gca,'xlim',[1 50])
set(gca,'ylim',[0 3.5])
box off

subplot(426)
plot(KaecEmStd)
title('std')
set(gca,'xlim',[1 50])
set(gca,'ylim',[0 YLIM_FOR_STD])
box off




subplot(427)
plot(KdiffMean)
%errorbar(KdiffMean,KdiffStd)
title('mean of Klin and Kref')
set(gca,'xlim',[1 50])
set(gca,'ylim',[0 3.5])
box off
subplot(428)
plot(KdiffStd)
title('std')
set(gca,'xlim',[1 50])
set(gca,'ylim',[0 YLIM_FOR_STD])
box off

subplot(421)
legend (sprintf('\\sigma_{in} = %.2f',Rs(1)),...
        sprintf('\\sigma_{in} = %.2f',Rs(2)),...
        sprintf('\\sigma_{in} = %.2f',Rs(3)),...
        sprintf('\\sigma_{in} = %.2f',Rs(4)))
legend boxoff
    

set(gcf, 'paperposition', [0 0 15, 12]) 
set(gcf, 'papersize', [15 12]) 
filename = sprintf('%s_linear_kernels_EM_vs_AEC_100.pdf',DATA)
saveas(1,filename)

%% replot for each input level
clf;
M=EM_FIT.EM.M;
for i=1:length(Rs)
    subplot(2,2,i); hold on
    %plot([KemMean(:,i) KaecMean(:,i) KdiffMean(1:50,i)])
    errorbar(KaecMean(:,i),KaecStd(:,i),'g')
    errorbar(KemMean(:,i),KemStd(:,i))
    errorbar(KaecEmMean(:,i),KaecEmStd(:,i),'r')
    
    errorbar(KdiffMean(1:M,i),KdiffStd(1:M,i),'k')
    
    title (sprintf('Ke estimate with \\sigma_{in}=%.2f',Rs(i)))
    
    set(gca,'xlim',[1 M])
    set(gca,'ylim',[0 3.5])

    box off
    
    if i==1
        legend ('AEC','EM','AEC-EM','ref');
        legend boxoff
    end
end

set(gcf, 'paperposition', [0 0 15, 12]) 
set(gcf, 'papersize', [15 12]) 
filename = sprintf('%s_linear_kernels_EM_vs_AEC_each_input100.pdf',DATA)
saveas(1,filename)
    
    
%% Question :is there difference in estimate as iteration inreases?

clear

AEC_EM_FIT0=load('../2012-07-08_AEC_DEC_itr50/fit_CellBsub_rep1-10_N5100_itr50_trial32.mat')
KaecEmMean0 = reshape(mean(AEC_EM_FIT0.KeEsts,2),[],size(AEC_EM_FIT0.KeEsts,3));
KaecEmStd0 = reshape(std(AEC_EM_FIT0.KeEsts,[],2),[],size(AEC_EM_FIT0.KeEsts,3));

AEC_EM_FIT=load('../2012-07-08_AEC_DEC_itr100/fit_CellBsub_rep1-10_N5100_itr100_trial32.mat')
KaecEmMean = reshape(mean(AEC_EM_FIT.KeEsts,2),[],size(AEC_EM_FIT.KeEsts,3));
KaecEmStd = reshape(std(AEC_EM_FIT.KeEsts,[],2),[],size(AEC_EM_FIT.KeEsts,3));


AEC_EM_FIT2=load('../2012-07-08_AEC_DEC_itr1000/fit_CellBsub_rep1-10_N5100_itr1000_trial32.mat')
KaecEmMean2 = reshape(mean(AEC_EM_FIT2.KeEsts,2),[],size(AEC_EM_FIT2.KeEsts,3));
KaecEmStd2 = reshape(std(AEC_EM_FIT2.KeEsts,[],2),[],size(AEC_EM_FIT2.KeEsts,3));



clf;
Rs=AEC_EM_FIT.Rs;
M=AEC_EM_FIT.EM.M;
for i=1:length(Rs)
    subplot(2,2,i); hold on
    errorbar(KaecEmMean0(:,i),KaecEmStd0(:,i),'k')
    errorbar(KaecEmMean(:,i),KaecEmStd(:,i),'b')
    errorbar(KaecEmMean2(:,i),KaecEmStd2(:,i),'r')
    
%     plot(KaecEmMean0(:,i),'k')
%     plot(KaecEmMean(:,i),'b')
%     plot(KaecEmMean2(:,i),'r')
    
%     plot(KaecEmStd0(:,i),'k')
%     plot(KaecEmStd(:,i),'b')
%     plot(KaecEmStd2(:,i),'r')
    
    %errorbar(KdiffMean(1:M,i),KdiffStd(1:M,i),'k')
    
    title (sprintf('Ke estimate with \\sigma_{in}=%.2f',Rs(i)))
    
    set(gca,'xlim',[1 M])
%     set(gca,'ylim',[0 3.5])

    box off
    
    if i==1
        legend ('50 itr', '100 itr','1000 itr');
        legend boxoff
    end
end

set(gcf, 'paperposition', [0 0 15, 12]) 
set(gcf, 'papersize', [15 12]) 
filename = sprintf('%s_linear_kernels_AEC-EM_itr_100_vs_1000.pdf', AEC_EM_FIT.EXP.DATA)
saveas(1,filename)
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% combine results by different method and plot in a same plot (June 10, 2012)
% load individual results and plot by running the above block

clear
load ../2012-06-09_multi_tap_GD/fit_HHspike_GD_rep1-50_dim1_N5100_itr50_trial32.mat Rs GD paramEstTrueV paramEstLinRegMulti mseSingleTap
ESTIMATE_MULTI_TAP_GD = 1

load ../2012-06-03_multi_tap_with_lin_reg/fit_HHspike_rep1-50_dim1_N5100_itr50_trial32.mat
% AtEsts paramEstMultis
ESTIMATE_MULTI_TAP_EM = 1

%%
clear
load ../2012-06-09_multi_tap_GD/fit_HHspike_GD_rep1-50_dim2_N5100_itr50_trial32.mat Rs GD paramEstTrueV paramEstLinRegMulti mseSingleTap
ESTIMATE_MULTI_TAP_GD = 1

load ../2012-06-03_multi_tap_with_lin_reg/fit_HHspike_rep1-50_dim2_N5100_itr50_trial32.mat
% AtEsts paramEstMultis
ESTIMATE_MULTI_TAP_EM = 1

%%
clear
load ../2012-06-09_multi_tap_GD/fit_HHspike_GD_rep1-50_dim3_N5100_itr50_trial32.mat Rs GD paramEstTrueV paramEstLinRegMulti mseSingleTap
ESTIMATE_MULTI_TAP_GD = 1

load ../2012-06-03_multi_tap_with_lin_reg/fit_HHspike_rep1-50_dim3_N5100_itr50_trial32.mat
% AtEsts paramEstMultis
ESTIMATE_MULTI_TAP_EM = 1

%%
clear
load ../2012-06-09_multi_tap_GD/fit_HHspike_GD_rep1-50_dim4_N5100_itr50_trial32.mat Rs GD paramEstTrueV paramEstLinRegMulti mseSingleTap
ESTIMATE_MULTI_TAP_GD = 1

load ../2012-06-03_multi_tap_with_lin_reg/fit_HHspike_rep1-50_dim4_N5100_itr50_trial32.mat
% AtEsts paramEstMultis
ESTIMATE_MULTI_TAP_EM = 1

%%
clear
load ../2012-06-09_multi_tap_GD/fit_HHspike_GD_rep1-50_dim5_N5100_itr50_trial32.mat Rs GD paramEstTrueV paramEstLinRegMulti mseSingleTap
ESTIMATE_MULTI_TAP_GD = 1

load ../2012-06-03_multi_tap_with_lin_reg/fit_HHspike_rep1-50_dim5_N5100_itr50_trial32.mat
% AtEsts paramEstMultis
ESTIMATE_MULTI_TAP_EM = 1















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% To compare roots for different dimensions  (June 10, 2012)


clear

%% plot roots 
% unit circle;
clear
%load('../2012-06-03_multi_tap_with_lin_reg/fit_HHspike_rep1-50_dim1_N5100_itr50_trial32.mat','Rs')
load ../2012-07-09_CellB_AEC_DEC/fit_CellBspike_rep1-10_dim1_N5100_itr100_trial32.mat Rs EXP


theta = 0:.01:(2*pi);
UC = [cos(theta); sin(theta)];


dims=1:4
numD = length(dims);
numRow=numD
numCol=length(Rs);
clf;

log_mseEM_mean = zeros(numD,length(Rs));
log_mseEM_std = zeros(numD,length(Rs));
qHatEM_mean = zeros(numD,length(Rs));
qHatEM_std = zeros(numD,length(Rs));


for dim = dims
    %filename = sprintf('2012-05-17_HH_nke4_additional/fit_HHspike_rep1-50_dim%d_N5100_itr500_trial32.mat',dim);
    %filename = sprintf('fit_HHspike_rep1-50_dim%d_N5100_itr50_trial32.mat',dim)
    
    
    % load EM results 
    filename = sprintf('../2012-07-09_CellB_AEC_DEC/fit_CellBspike_rep1-10_dim%d_N5100_itr100_trial32.mat',dim)
    load(filename, 'AtEsts', 'dimX','mseMultiTap','paramEstMultis')
    
%     % load GD results 
%     filename = sprintf('../2012-07-09_CellB_AEC_DEC/fit_CellBspike_rep1-10_dim%d_N5100_itr100_trial32.mat',dim)
%     load(filename, 'GD');   %'AtEsts', 'dimX','mseMultiTap','paramEstMultis')
    
    
% %     % load EM results 
% %     filename = sprintf('../2012-06-03_multi_tap_with_lin_reg/fit_HHspike_rep1-50_dim%d_N5100_itr50_trial32.mat',dim)
% %     load(filename, 'AtEsts', 'dimX','mseMultiTap','paramEstMultis')
% %     
% %     % load GD results 
% %     filename = sprintf('../2012-06-09_multi_tap_GD/fit_HHspike_GD_rep1-50_dim%d_N5100_itr50_trial32.mat',dim)
% %     load(filename, 'GD');   %'AtEsts', 'dimX','mseMultiTap','paramEstMultis')

    for i=1:length(Rs)
        R = Rs(i);
        subplot(numRow,numCol,i+(dim-1)*numCol); hold on
        %for rep = 1:size(AtEsts,2);
        for rep = 1:1
            % plot root of EM
            root = roots([1;-AtEsts(1:dimX,rep,i)]);
            plot(real(root),imag(root), 'ob')
            
            % plot root of GD
            if exist('GD','var')
                root = roots([1;-GD.AtEst(1:dimX,rep,i)]);
                plot(real(root),imag(root), 'sg')
            end
            
        end
        % plot unit circle
        plot (UC(1,:), UC(2,:), ':k')
        box off
        
        
        
        
        if dim==1
            %title (sprintf('R=%.1f',R))
            title (sprintf('\\sigma_{in}=%.2f',R))
        end
        
        if i==1
            ylabel (sprintf('tap=%d',dim))
        end
        
        if dim==1 & i==1 & exist('GD','var')
            legend ('EM', 'GD', 'location','NW'); legend boxoff
        end
    end
    
    % MSE
    log_mseEM_mean(dim,:) = mean(log2(mseMultiTap));
    log_mseEM_std(dim,:) = std(log2(mseMultiTap));

    % check q
    qHatEM_mean(dim,:) = mean(paramEstMultis(end-1,:,:),2);
    qHatEM_std(dim,:) = std(paramEstMultis(end-1,:,:),[],2);
    
    
end
%

set(gcf, 'paperposition', [0 0 numRow*2 numCol*1.5]) 
set(gcf, 'papersize', [numRow*2 numCol*1.5]) 
filename = sprintf('%s_roots.pdf',EXP.baseFilename)
saveas(1,filename)
        








%% plot MSE and Q  (??)
clf
subplot(121)
errorbar(repmat(log2(Rs),numD,1)',qHatEM_mean',qHatEM_std','o-'); xlabel('log2(R)')
if exist('q','var')
    hold on
%    plot(log2(Rs),repmat(q,1,length(Rs)),'k--+')
end
legend('tap=1','tap=2','tap=3','tap=4','tap=5','true','location','nw');legend boxoff
title ('Qhat'); box off
subplot(122)
errorbar(repmat(log2(Rs),numD,1)',log_mseEM_mean',log_mseEM_std','o-'); xlabel('log2(R)')
hold on; plot(log2([Rs(1) Rs(end)]),log2([Rs(1) Rs(end)]),'k+--')
%axis equal
title ('log2(MSE)'); box off



set(gcf, 'paperposition', [0 0 8, 3]) 
set(gcf, 'papersize', [8 3]) 
filename = sprintf('%s_q_mse.pdf',EXP.baseFilename)
saveas(1,filename)




%% 




























load mseEMKeLin_N5100.mat mseEMKeLin
LINE_TYPE = '--x';
hE = errorbar(Rs(1:numR), mean(mseEMKeLin), std(mseEMKeLin),LINE_TYPE,'lineWidth',LINE_WIDTH); hold on
        set(gca,'Xscale','log')
        set(gca,'Yscale','log')

        % adjust error bar width
        hE_c          = get(hE     , 'Children'    );
        errorbarXData = get(hE_c(2), 'XData'       );
        BAR_WIDTH_INC =  1.2;
        errorbarXData          = ...
        get(hE_c(2), 'XData'       );
        errorbarXData(4:9:end) = errorbarXData(1:9:end) / BAR_WIDTH_INC;
        errorbarXData(7:9:end) = errorbarXData(1:9:end) / BAR_WIDTH_INC;
        errorbarXData(5:9:end) = errorbarXData(1:9:end) * BAR_WIDTH_INC;
        errorbarXData(8:9:end) = errorbarXData(1:9:end) * BAR_WIDTH_INC;
        set(hE_c(2), 'XData', errorbarXData);

        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot multiple Ns on top of each other 
%numR = length(Rs);
LINE_WIDTH = 1;
%load fit_HHsub_rep1-10_N1100_itr1000_trial16.mat
%load fit_HHsub_rep1-10_N5100_itr2000_trial32
%load fit_HHsub_rep1-10_N5100_itr10000_trial32
%load fit_HHsub_rep1-50_N5100_itr10000_trial32
%load fit_HHsub_rep1-10_N5100_itr2000_trial32
load fit_HHsub_rep1-50_N5100_itr1000000_trial32
%load fit_HHspike_rep1-10_N5100_itr5000_trial32  % for spiking data
LINE_TYPE = 'k'
LINE_TYPE2 = 'k--x'
PLOT_REFERENCE = 0
clf
% % %%
% % %load fit_HHsub_rep1-10_N6000_itr1000_trial16.mat
% % %load fit_HHsub_rep1-10_N20100_itr2000_trial32
% % LINE_TYPE = 'b'
% % PLOT_REFERENCE = 1
%%
%load fit_HHsub_rep1-10_N11000_itr1000_trial16.mat
%load fit_HHsub_rep1-10_N10100_itr2000_trial32
%load fit_HHsub_rep1-10_N10100_itr10000_trial32
%load fit_HHsub_rep1-50_N10100_itr10000_trial32
%load fit_HHsub_rep1-10_N10100_itr2000_trial32
load fit_HHsub_rep1-50_N10100_itr1000000_trial32
%load fit_HHspike_rep1-10_N10100_itr5000_trial32  % for spiking data


LINE_TYPE = 'b'
LINE_TYPE2 = 'b--x'
PLOT_REFERENCE = 0
%%
%load fit_HHsub_rep1-10_N21000_itr1000_trial16.mat
%load fit_HHsub_rep1-10_N20100_itr2000_trial32
%load fit_HHsub_rep1-10_N20100_itr10000_trial32
load fit_HHsub_rep1-50_N20100_itr1000000_trial32

LINE_TYPE = 'g'
LINE_TYPE2 = 'g--x'
PLOT_REFERENCE = 1

%% 
%legend ('N=1000','N=5000','N=10000','N=20000'); legend boxoff
% legend ('N=5000','N=10000','N=20000','location','SE'); legend boxoff

subplot(121)
legend ('N=5000','N=10000','N=20000','location','NW'); legend boxoff
subplot(122)
legend ('N=5000','N=5000 w/o sm', 'N=10000', 'N=10000 w/o sm','N=20000', 'N=20000 w/o sm','location','NW'); legend boxoff
%legend ('N=5000','N=20000'); legend boxoff

%% save 
if ~exist('mseMultiTap','var')
%     subplot(121);set(gca,'xlim', 10.^[-2.1 0.1])
%     subplot(122);set(gca,'xlim', 10.^[-2.1 0.1])
    set(gcf, 'paperposition', [0 0 9 4]) 
    set(gcf, 'papersize', [9 4]) 
    %legend ('N=5000','N=10000','N=20000','Ke w/o smooth','location','SE'); legend boxoff
else
    set(gcf, 'paperposition', [0 0 5 4]) 
    set(gcf, 'papersize', [5 4]) 
    legend ('N=5000','N=10000','N=20000','Ke w/o smooth','location','NW'); legend boxoff
end
    
filename = sprintf('%s%s_rep1-%d_MSE_Ns_itr%d_trial%d.pdf',EXP.DATA,EXP.METHOD,EXP.numRep,EM.max_iter,EM.num_trial)    
saveas(1, filename)


%% load AEC
%fit_HHAEC_rep1-10_N5100_itr10000_trial32
clear
clf
load fit_HHsub_rep1-10_N5100_itr10000_trial32
LINE_TYPE ='b'

%%
load fit_HHsub_rep1-10_N10100_itr10000_trial32
LINE_TYPE ='g'
%%
load fit_HHAEC_rep1-10_N5100_itr10000_trial32
LINE_TYPE ='k'

%%
clf

%%
subplot(121); hold on
%plot(Rs,mseKe')
errorbar(Rs,mean(mseKe),std(mseKe),LINE_TYPE)
set(gca,'xscale','log')
set(gca,'yscale','log')

subplot(122); hold on
%plot(Rs,mse')
errorbar(Rs,mean(mse),std(mse),LINE_TYPE)
set(gca,'xscale','log')
set(gca,'yscale','log')



%%

legend ('EM N=5000', 'EM N=10,000', 'AEC N=5000')
saveas(1,'comparision.pdf')

%% save N=5000 vs N=10,000
subplot(121)
legend ('N=50000','N=10000'); legend boxoff
subplot(122)
legend ('smoothing','w/o smoothing','location','NW'); legend boxoff
%legend ('smoothing (N=5000)','w/o smoothing (N=5000)','smoothing (N=10000)','w/o smoothing (N=10000)'); legend boxoff
saveas(1,'comparision for diff Ns.pdf')

%% 

clear
clf
load fit_HHsub_rep1-50_N5100_itr1000000_trial32
LINE_TYPE ='b'

%%
load fit_HHsub_rep1-50_N10100_itr100000_trial32
LINE_TYPE ='g'
%%
load fit_HHAEC_rep1-10_N5100_itr10000_trial32
LINE_TYPE ='k'



%%
legend ('EM N=5000', 'EM N=10,000', 'AEC N=5000')
saveas(1,'comparision_rep50.pdf'); legend boxoff









%% 
set(hE1                            , ...
  'LineStyle'       , 'none'      , ...
  'Marker'          , '.'         , ...
  'Color'           , [.3 .3 .3]  );







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's plot estimated params (Apr. 26, 2012)


%% fit spiking data with single tap
clear
load fit_HHspike-single_rep1-10_N5100_itr5000_trial32.mat
%load fit_HHspike-single_rep1-50_N5100_itr5000_trial32.mat
whos

%% fit spiking data with multi tap 
clear
load fit_HHspike_rep1-10_dim2_N5100_itr5000_trial32
%load fit_HHspike_rep1-10_dim3_N5100_itr5000_trial32
whos

% paramTrue is not correct ..

%%
switch EXP.METHOD
    case {'sub','spike-single'}
        numR = length(Rs);
        mu = mean(paramEsts,2);
        mu = reshape (mu, [], numR);
        sig = std(paramEsts,0,2);
        sig = reshape (sig, [], numR);

        muParam = mu(1:3,:);
        muQR = mu(end-1:end,:);
        sigParam = sig(1:3,:);
        sigQR = sig(end-1:end,:);
        
        

        muKe = reshape(mean(KeEsts,2),[], length(Rs));
        sigKe = reshape(std(KeEsts,0,2),[],length(Rs));
        
        numCol=3;
        
        X_TICK_LABEL = {'\alpha','\beta', 'V_{rev}'};
        
    case 'spike'
        dimX = size(AtEsts,1)
        switch dimX
            case 2
                X_TICK_LABEL = {'alpha1','alpha2'}
            case 3
                X_TICK_LABEL = {'alpha1','alpha2','alpha3'}
                    
        end
        
        numR = length(Rs);
        % calc mean
        muA = mean(AtEsts,2);
        muA = reshape (muA, [], numR);
        muEst = mean(paramEstMultis,2);
        muEst = reshape (muEst, [], numR);        
        % combine
        muParam = [muA; muEst(1:end-2,:)];
        muQR = muEst(end-1:end,:);
        muQR = muEst(end-1:end,:);
        
        % calc std
        sigA = std(AtEsts,0,2);
        sigA = reshape (sigA, [], numR);
        sigEst = std(paramEstMultis,0,2);
        sigEst = reshape (sigEst, [], numR);
        % combine
        sigParam = [sigA; sigEst(1:end-2,:)];
        sigQR = sigEst(end-1:end,:)
        
        numCol=4;
end


%
clf

for idxR = 1:numR
    
    %
    subplot(numR,numCol,numCol*(idxR-1)+1)
%     idx=1:3;
    errorbar(muA(:,idxR), sigA(:,idxR),'x:')
    switch EXP.DATA
        case 'HH'   % have only true Vrev and R
%             hold on
%             plot([NaN NaN paramTrue.Vo],'r+--')     % doesn't have much meaning
%     if exist (EXP.TRUE)
%         hold on;
%         plot([paramTrue.alpha paramTrue.beta paramTrue.gamma], 'r+--')
%     end
    end
%             
    set(gca,'xtick',1:length(X_TICK_LABEL))
    set(gca,'xticklabel',X_TICK_LABEL)
    %title (sprintf('\\alpha , \\beta, V_{rev} for R=%.2f',Rs(idxR)))
    title (sprintf('Estimate for R=%.2f',Rs(idxR)))
    box off
    
    
    % plot beta separately 
     subplot(numR,numCol,numCol*(idxR-1)+2)
%     idx=1:3;
    errorbar(muParam(end-1,idxR), sigParam(end-1,idxR),'x')
    switch EXP.DATA
        case 'HH'   % have only true Vrev and R
%             hold on
%             plot([NaN NaN paramTrue.Vo],'r+--')     % doesn't have much meaning
%     if exist (EXP.TRUE)
%         hold on;
%         plot([paramTrue.alpha paramTrue.beta paramTrue.gamma], 'r+--')
%     end
    end
    set(gca,'xlim', [0.8 1.2])
%             
    set(gca,'xtick',1)
    set(gca,'xticklabel','beta')
    box off
    
    % plot Vrev separately 
     subplot(numR,numCol,numCol*(idxR-1)+3)
%     idx=1:3;
    errorbar(muParam(end,idxR), sigParam(end,idxR),'x')
    switch EXP.DATA
        case 'HH'   % have only true Vrev and R
%             hold on
%             plot([NaN NaN paramTrue.Vo],'r+--')     % doesn't have much meaning
%     if exist (EXP.TRUE)
%         hold on;
%         plot([paramTrue.alpha paramTrue.beta paramTrue.gamma], 'r+--')
%     end
    end
    set(gca,'xlim', [0.8 1.2])
%             
    set(gca,'xtick',1)
    set(gca,'xticklabel','Vo')
    box off
    
    
    %
    subplot(numR,numCol,numCol*(idxR-1)+4)
    errorbar(muQR(:,idxR), sigQR(:,idxR));
    if (EXP.TRUE)
        hold on;
        %plot([paramTrue.Q paramTrue.R], 'r+--')
        plot([Qs Rs(idxR)], 'r+:')
        
    end    
    set(gca,'xtick',1:2)
    set(gca,'xticklabel',{'Q','R'})
%     title ('Q and R')
    box off
    if idxR==1
        legend ('Estimate', 'True','location','NE'); legend boxoff
    end

    
    %
    %
    switch EXP.METHOD
        case {'sub','spike-single'}     % ke is estimated for those cases
        subplot(numR,numCol,numCol*(idxR-1)+3)
        idx=4:5;
        errorbar(muKe(:,idxR),sigKe(:,idxR)); 
        if (EXP.TRUE)
            hold on
            plot (.5*[.8 7 3 .3],'r+:')
        end
        title ('K_e')
        box off
    end
    
    %
%     subplot(2,2,4)
%     idx=EM.M+6:EM.M+7;
%     errorbar(mu(idx,idxR), sig(idx,idxR));
%     if (EXP.TRUE)
%         hold on;
%         %plot([paramTrue.Xo paramTrue.P], 'r+--')       % FIXME
%         plot ([-60 1], 'r+--')
%     end    
%     title ('X_o and P_o')


end

set(gcf, 'paperposition', [0 0 12, 3*length(Rs)]) 
set(gcf, 'papersize', [12 3*length(Rs)]) 
switch EXP.METHOD
    case {'sub','spike-single'}
            
        filename = sprintf('param_%s_%s_rep1_%d_N%d_itr%d_trial%d.pdf',EXP.DATA, EXP.METHOD,EXP.numRep, EM.T,EM.max_iter,EM.num_trial)
    case {'spike'}
        filename = sprintf('param_%s_%s_rep1_%d_N%d_dim%d_itr%d_trial%d.pdf',EXP.DATA, EXP.METHOD,  EXP.numRep,  EM.T, dimX, EM.max_iter,EM.num_trial)
end
saveas(1, filename)

%%
% paramEst.Cm = dt/paramEst.beta/sum(paramEst.Ke);
% paramEst.gL = (1-paramEst.alpha)/paramEst.beta/sum(paramEst.Ke);
% paramEst.Vo = paramEst.gamma ./ (1-paramEst.alpha);

[numParam numRep numRs] = size(paramEsts);

% accuracy of parameters 
% calc Cm
%dt = .1;  %M=7;
Cm = dt./paramEsts(2,:,:)./sum(paramEsts(4:3+EM.M,:,:),1);
Cm = reshape(Cm, numRep, numRs);
% calc gL
gL = (1-paramEsts(1,:,:))./paramEsts(2,:,:)./sum(paramEsts(4:3+EM.M,:,:),1);
gL = reshape(gL, numRep, numRs);
% calc Vo
Vo = paramEsts(3,:,:) ./ (1-paramEsts(1,:,:));
Vo = reshape(Vo, numRep, numRs);
%%
clf
idxToShow = 1:numRs;
subplot(311)
errorbar(Rs(idxToShow), mean(Cm(:,idxToShow)), std(Cm(:,idxToShow)));
if (EXP.TRUE)
    hold on; semilogx (Rs(idxToShow), repmat(1,1,length(idxToShow)),'r--x')
    set(gca,'ylim', [0 2])
    legend ('estimate', 'true')
end
title ('C_m')

subplot(312)
errorbar(Rs(idxToShow), mean(gL(:,idxToShow)), std(gL(:,idxToShow)))
if (EXP.TRUE)
    hold on; semilogx (Rs(idxToShow), repmat(0.3,1,length(idxToShow)),'r--x')
    set(gca,'ylim', [0 0.5])
    legend ('estimate', 'true')
end
title ('g_L')

subplot(313)
errorbar(Rs(idxToShow), mean(Vo(:,idxToShow)), std(Vo(:,idxToShow)))
if (EXP.TRUE)
    hold on; semilogx (Rs(idxToShow), repmat(10.6,1,length(idxToShow)),'r--x')
    set(gca,'ylim', [-30 50])
    legend ('estimate', 'true')
end
title ('V_o')


filename = sprintf('paramEsts_%s_R%.2f_Mstep%d%d_T%d.mat',EXP.DATA, R, EM.MstepMethod, EM.MstepConstraint, EM.T)
save (filename, 'Cm', 'gL', 'Vo', 'Rs');
%save paramEsts Cm gL Vo Rs

filename = sprintf('paramEsts_%s_R%.2f_Mstep%d%d_T%d.pdf',EXP.DATA, R, EM.MstepMethod, EM.MstepConstraint, EM.T)
saveas (1, filename)



return

















%% load and plot results
% % for i=1:length(Rs)
% %     
% %     filename = sprintf('fit_%s_R%.2f_Mstep%d_T%d.mat',EXP.METHOD, R, EM.MstepMethod, EM.T)
% %     save (filename, 'Rs', 'paramEsts', 'params', 'mse')
% % end    
    
%load fit_sub_R1.00_Mstep2_T5000.mat
%load ../fit_EM_Mstep0/fit_sub_R1.00_Mstep2_T5000.mat
load ../fit_EM_Mstep0/fit_sub_R10.00_Mstep0_T5000.mat
numR = length(Rs)
clf;
if (size(mse,1)~=1)
    errorbar(Rs(1:numR), mean(mse(:,1:numR)), std(mse(:,1:numR)))
else
    semilogx(Rs, mse,'+-')
end
ylabel('mse')
xlabel('R')

hold on
load ../fit_EM_Mstep2_10/fit_sub_R1.00_Mstep2_T5000.mat
numR = 5
if (size(mse,1)~=1)
    errorbar(Rs(1:numR), mean(mse(:,1:numR)), std(mse(:,1:numR)),'r')
else
    semilogx(Rs, mse,'+-',r')
end
ylabel('mse')
xlabel('R')

legend ('Mstep0', 'fmincon')


%filename = sprintf('MSE_%s_Mstep%d_T%d.pdf',EXP.METHOD,EM.MstepMethod,EM.T)
%saveas(1, filename)



    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of routine here
% below are the quick code done manually to check results. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% compare with AEC (Jan. 08, 2010)

%load ../fit_AEC_hh/aec_MSE_Rs
load ../fit_AEC_hh_T5000/aec_MSE_Rs MSE
load fit_HHsub_T5000 mse Rs
clf
semilogx(Rs, [MSE mse'],'x-')
legend ('AEC', 'EM');
ylabel('MSE')
xlabel('R')
filename = sprintf('MSE_AEC v.s. EM.pdf')
saveas(1,filename)

%% Compare AEC, EM, EM+denoise (Jan. 09, 2010)
% compensate V by Y-K*II where K is given by EM
%load fit_HHsub_T5000 mse Rs
clear EMK
EMK.MODE='sub';
EMK.T = 5000;
EMK.mse = [];
for R=Rs;
    % load data generated by HH
    filename = sprintf('../Data/HH%s_R%.2f.mat',EMK.MODE,R);
    load (filename)
    
    % clip data length to be AEC.T
    Iapp = Iapp(1:EMK.T);
    II = II(:,1:EMK.T);
    V = V(1:EMK.T);  
    X = V';
    Ve = Ve(1:EMK.T);
    Y = Y(1:EMK.T);  
    
    % load estimated kernel 
    filename = sprintf('fit_HH%s_R%.2f.mat',EMK.MODE,R)
    load (filename, 'paramEst', 'mse')
    
    % compensate 
    EMK.Ve = paramEst.Ke*II;
    EMK.Vm = Y-EMK.Ve;
    EMK.mseVe = (Ve-EMK.Ve)*(Ve-EMK.Ve)'/EMK.T
    EMK.mseVm = (X-EMK.Vm)*(X-EMK.Vm)'/EMK.T
    
    EMK.mse = [EMK.mse; EMK.mseVm];    
end

filename = sprintf('EMK_HH%s_T%d.mat',EMK.MODE,EMK.T)
save (filename, 'Rs','EMK')

clf
semilogx(Rs, EMK.mse,'x-')
ylabel('MSE')
xlabel('R')
filename = sprintf('EMK_HH%s_T%d.pdf',EMK.MODE,EMK.T)
saveas(1,filename)

%% draw and compare all three cases  (Jan. 15, 2010)
load ../fit_AEC_hh_T5000/aec_MSE_Rs MSE
load fit_HHsub_T5000 mse Rs
load EMK_HHsub_T5000 EMK
clf
semilogx(Rs, [MSE EMK.mse mse'],'x-','linewidth',3,'MarkerSize',12)
legend ('AEC', 'EM', 'EM+denoise');
ylabel('MSE')
%xlabel('R')
xlabel('Measurement noise (\sigma^2_{rec})')
title('Error in estimating V_m')
%filename = sprintf('MSE_AEC v.s. EMs.pdf')
%saveas(1,filename)

set(gcf, 'paperposition', [.25 .25 4 3])    % save with thicker line and larger text
print -dpng MSE_comparison.pdf


%% check estimated cours
idx = 4;
clf
subplot(311)
plot (Iapps(1:T,idx),'r');hold on
subplot(312)
plot (Ys(1:T,idx))

title ('measured')
subplot(313)
plot (Vs(1:T,idx),'r'); hold on
plot(paramEst(1,idx).Xs);
legend ('True', 'Estim')






%% Plot time courses (Jan. 15, 2010)
idx = 1:1000;
clf;
subplot(311);
plot(dt*idx,Iapp(idx),'r');hold on
plot(dt*idx,Iinj(idx));
legend('I_{app}','I_{inj}');
xlabel('ms')
ylabel('uA/cm^2')


subplot(312);
plot(tt(idx),V(idx));
title('V_m');
xlabel('ms')
ylabel('mV')

subplot(313);
plot(tt(idx),[m(idx) n(idx) h(idx)]);
%plot(tt,[m n h]);
title('Gating variables');
xlabel('ms');
ylabel('1/ms')
legend ('m', 'n', 'h')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate alpha,beta,Ke,Q,R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EM.num_trial = 10;
EM.max_iter = 100;

% Initial guess of kernel using simple linear regression
Klin = Y'\II';
Klin = flipud(Klin);

clf; 
clear paramInit paramEst paramEstim LLtrace;
LLend = [];
XXs = zeros(EM.num_trial,T);
for trial = 1:EM.num_trial
    disp(sprintf('%dth trial',trial))
    
    % initialize
    paramInit(trial).alpha = rand();
    paramInit(trial).beta = 0.02*rand();
    paramInit(trial).gamma = -.65;
    paramInit(trial).C = 1;
    paramInit(trial).Ke = Klin + randn(1,M);%rand(size(Ke));
    paramInit(trial).Q = rand();
    paramInit(trial).R = 2*R*rand();
    paramInit(trial).Xo = randn();
    paramInit(trial).Po = 2*rand();
    
%     % start some of the value from the original value (for debug)
%     paramInit(trial).alpha = alpha;
%     paramInit(trial).beta = beta;
%     paramInit(trial).Ke = Ke;
%     paramInit(trial).Q = Q;
%     paramInit(trial).R = R;
%     paramInit(trial).Xo = Xo;
%     paramInit(trial).Po = Po;
       
    % estimate parameters
    %[paramEstim(trial) LL, XXs(trial,:)] = em_kalman_abc(Y, II, paramInit(trial), EM.max_iter);
    [paramEstim(trial) LL, XXs(trial,:)] = em_kalman_nn(Y, II, paramInit(trial), EM.max_iter);     % with non-negative constraint on Ke
    LLtrace(trial).LL = LL;   
    LLend = [LLend LL(end)];
    
    disp(sprintf('  LL=%.2e', LL(end)))
    
    subplot(521); hold on;
    semilogy(real(LL));
    %title ('Log likelihood')
    title ('L(\Theta)')
    xlabel ('Iteration')
    subplot(522); hold on;
    plot(trial, LL(end),'x');
    title ('Final Log likelihood')
    xlabel ('Trial')
    drawnow;
end
% choose parameter set with maximum LL
[LLmax idx] = max(LLend);
paramEst = paramEstim(idx);
paramEst.Xs = XXs(idx,:);


Ifil = filter(Ke, 1, Iinj);
IfilEst = filter(paramEst.Ke, 1, Iinj);



paramEst.Cm = dt/paramEst.beta/sum(paramEst.Ke);
paramEst.gL = (1-paramEst.alpha)/paramEst.beta/sum(paramEst.Ke);
paramEst.Vo = paramEst.gamma ./ (1-paramEst.alpha);


% original parameters
gL = .1;
Cm =1;
Vo = -65;

alpha = 1-dt*gL/Cm;
beta = gL/Cm/sum(Ke);
gamma = gL*dt*Vo/Cm;


% save result
filename = sprintf('estimation_nn_R%.0f_itr%d.mat', R, EM.max_iter);
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
filename = sprintf('estimation_nn_R%.0f_itr%d.pdf',R,EM.max_iter);
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

filename = sprintf('generated_Q%.0f_R%.0f_T%d_itr%d.pdf', Q,R,T,EM.max_iter);
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
VeAEC = AEC.Keff(1:EM.M)'*II;
VmAEC = Y-VeAEC;
mseVeAEC = (Ve-VeAEC)*(Ve-VeAEC)'/T
mseVmAEC = (X-VmAEC)*(X-VmAEC)'/T

% 2) compensage by Ke estimated by EM
VeEMKe = paramEst.Ke*II;
VmEMKe = Y-VeEMKe;
mseVeEMKe = (Ve-VeEMKe)*(Ve-VeEMKe)'/T
mseVmEMKe = (X-VmEMKe)*(X-VmEMKe)'/T

% 3) smoothing by alpha, beta, gamma and Ke estimated by EM
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
legend ('True', 'AEC', 'EM')
xlabel('ms')
ylabel('mV')
title (sprintf('Electrode voltage V_e (MSE_{AEC} =%.2f, MSE_{EM}=%.4f)',mseVeAEC, mseVeEMKe))
set(gca,'xlim',[0 100])
subplot(212);hold off
plot (dt*(1:T),X(1,:),'r');hold on
plot (dt*(1:T),VmAEC,'g');
plot (dt*(1:T),VmEMKe,'c--');
plot (dt*(1:T),paramEst.Xs(1,:));
legend ('True', 'AEC', 'EM Ke','EM smooth')
xlabel('ms')
ylabel('mV')
title (sprintf('Membrane voltage (MSE_{AEC}=%.2f, MSE_{EM Ke}=%.2f,MSE_{EM smooth}=%.2f)',mseVmAEC, mseVmEMKe, mseVmEMsm))
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

EM = paramEst
save est X Y II alpha beta C Ke Q R Xo Po EM

