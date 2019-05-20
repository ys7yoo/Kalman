clear
%filename = '../fit_LDS_AEC_Q1_1205/all_LDS_AEC_T10000.mat'
filename = '../fit_LDS_EM_Q.1_1205/all_LDS_EM_Mstep00_T10000_itr100.mat'
EM = load(filename)
filename = '../fit_LDS_AEC_Q.1_1205/all_LDS_AEC_T10000.mat'
AEC = load(filename)





%% 

Rs = EM.Rs;
Rs = 1;

%for repParam = 1:length(paramEsts)
for idxR = 1:length(Rs)
    
    % to read data
    load ../Data/LDS_Q0.01_R0.01.mat
    
    
    for repParam = 1
        paramAEC = AEC.paramEsts(repParam,idxR)
        paramEM = EM.paramEsts(repParam,idxR)
        
        A = paramEM.alpha;

        A = paramEM.alpha;
        B = [paramEM.beta*paramEM.Ke paramEM.gamma];
        C = 1;
        D = [paramEM.Ke 0];
        Q = paramEM.Q;
        R = paramEM.R;
        
        %for repTest = 1:length(paramEsts)
        M = length(paramEM.Ke);
        N = size(EM.Iapps,1);
        for repTest = 2
            Y = EM.Ys(:,repTest)';
            Vm = Vms(:,repTest)
            II=stackCols(EM.Iapps(:,repTest), M,0);
            UU = [II; ones(1,N)];

            Xo = Y(1);
            Po = max(squeeze(paramEM.Ps));
            %[Xp Pp Xf Pf Kf LL] = kalman_filt(Y, UU, A, B, C, D, Q, R, Xo, Po);

            
            %% correction by subtraction
            % 1) AEC
            Ve = paramAEC.Ke'*II;
            VmSubAEC = Y-Ve;
            mse_subAEC = (VmSubAEC-')*(VmSubAEC-Vm')'/length(VmSub)
            
            
            % 2) DEC
            Ve = paramEM.Ke*II;
            VmSub = Y-Ve;
            mse_sub = (VmSub-Vm')*(VmSub-Vm')'/length(VmSub)

                    
            
            %% now filter/smoother with the best paramEM
            [Xs Ps Pcs LL Xf Pf] = kalman_smth(Y, UU, A, B, C, D, Q, R, Xo, Po);
            mse_fltr = (Xf-Vm')*(Xf-Vm')'/length(Xf)
            mse_smth = (Xs-Vm')*(Xs-Vm')'/length(Xs)
            

            
            %% plot results

            clf;
            T = 500;
            subplot(211); hold on
            
            plot (dt*(1:T),VmSub(1:T),'k');
            plot (dt*(1:T),Xf(1:T),'b');
            plot (dt*(1:T),Xs(1:T),'g');
            plot (dt*(1:T),Vms(1:T,repTest),'r--');
            legend ('subtraction','filtering', 'smoothing', 'true')
            
            
%             subplot(312); hold on
%             plot (squeeze(Pf), 'b')
%             plot (squeeze(Ps), 'g')
            
            subplot(212); hold on
            plot (dt*(1:T), VmSub(1:T)-Vms(1:T,repTest)', 'k')
            plot (dt*(1:T), Xf(1:T)-Vms(1:T,repTest)', 'b')
            plot (dt*(1:T), Xs(1:T)-Vms(1:T,repTest)', 'g')
            title (sprintf('error (mse_{sub}=%.2f,mse_{fltr}=%.2f,mse_{smth}=%.2f)',mse_sub,mse_fltr,mse_smth))
            
            

        end
    end
end