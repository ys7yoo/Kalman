function Theta =m_step_opt(Y, U, Xs, Ps, Pcs, Xo, Po, Theta0, const)
% reimplement M-step with optimization function



if (nargin < 9)
    const = 0;
end
    

    

% function to be minimized
%lossDerivFunc=@(Theta)(calcNegLLdLL (Xs, Y, U, Ps, Pcs, Xo, Po, Theta));
lossDerivFunc=@(Theta)(calcNegLLparamLog(Xs, Y, U, Ps, Pcs, Xo, Po, Theta));



if (const == 0)
    % 1) fminunc
        
%     tic
%     % using gradient    => gradient is not working yet
%     options = optimset('display','on', 'MaxFunEvals', 1000, 'TolX', 1e-5, 'GradObj', 'on');
%     %options = optimset('display','off', 'TolFun', 1, 'MaxIter', 50, 'largescale', 'off', 'GradObj', 'on');
%     [Theta,negLL] = fminunc(lossDerivFunc,Theta0,options);
%     toc
    
    

    %tic;
    %options = optimset('GradObj', 'on', 'display','off', 'MaxFunEvals', 1000, 'TolX', 1e-5);
    options = optimset('GradObj', 'on', 'display','off', 'MaxFunEvals', 1000);
    %options = optimset('GradObj', 'off', 'display','off', 'MaxFunEvals', 1000);
    
    [Theta,negLL, ~, ~, g, H] = fminunc(lossDerivFunc,Theta0,options);
    %toc

   
% 
%     % using gradient    => not much better 
%     tic;
%     options = optimset('display','on', 'MaxFunEvals', 1000, 'TolX', 1e-5, 'GradObj', 'on');
%     %options = optimset('display','off', 'TolFun', 1, 'MaxIter', 50, 'largescale', 'off', 'GradObj', 'on');
%     [ThetaGrad,negLLGrad] = fminunc(lossDerivFunc,Theta0,options);
%     toc
%     
%     
%     [Theta0 Theta ThetaGrad]
%     [negLL negLLGrad]
    
    

    
    %options = optimset('display','off', 'TolFun', 1, 'MaxIter', 50, 'largescale', 'off', 'GradObj', 'on');
    %options = optimset('display','off', 'TolFun', 1, 'MaxIter', 50, 'largescale', 'off');   
else
    % 2) fmincon with non-negative constraint on Ke
    eps = 1e-5;
    Amax = 1e3;
    varmax = 1e3;
    kmin = -eps;
    kmax = 1e4;
    M = length(Theta0)-5;
    
    lb = [eps -eps -Amax kmin*ones(1,M)    -2*[1 1]]';
    ub = [1-eps Amax Amax kmax*ones(1,M)  Inf Inf]';
    
    % make sure initial Theta0 satisfies 
    Theta0 = max(Theta0,lb);
    Theta0 = min(Theta0,ub);
    
    
    %options = optimset('GradObj', 'on', 'display','on', 'MaxFunEvals', 1000, 'TolX', 1e-5);
    options = optimset('GradObj', 'on', 'display','on', 'MaxFunEvals', 1000);
    %options = optimset('GradObj', 'off', 'display','off', 'MaxFunEvals', 1000);
    
    
    %options = optimset('display', 'off', 'algorithm', 'active-set', 'maxfunevals', 10000);
%     options = optimset('display', 'iter', 'algorithm', 'active-set',  'GradObj', 'on');
    %options = optimset('display','off', 'TolFun', 1, 'MaxIter', 50, 'largescale', 'on');
    %options = optimset('display','off', 'TolFun', 1, 'MaxIter', 50, 'largescale', 'on', 'GradObj', 'on');    
    [Theta,negLL] = fmincon(lossDerivFunc,Theta0,[],[],[],[],lb,ub,[], options);


    
    %disp(sprintf('change in LL %f', LLstore-negLL));
    %keyboard;
    %LLstore = negLL;
end



DerivCheck(lossDerivFunc, Theta, options);
        