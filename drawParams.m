function drawParams(paramEstim, paramInit, paramTrue, lineColor)

if nargin <2
    paramInit = [];
end


if nargin <3
    paramTrue = [];
end

if nargin <4
    lineColor = [0 0 1];
end

if nargin <5
    highColor = [1 0 0];
end

HIGH_LIGHT_BEST_PARAM =1;


if isstruct(paramEstim)
    numParams = numel(paramEstim);
    nke = length(paramEstim(1).Ke);

% %     %convert to a matrix
% %     a = [paramEstim(:).alpha];
% %     b = [paramEstim(:).beta];
% %     v = [paramEstim(:).Vrev];
% %     ke = reshape([paramEstim(:).Ke],nke,[]);
% %     q = [paramEstim(:).Q];
% %     r = [paramEstim(:).R];
% %     
% %     params = [a;b;v;ke;q;r];

    % find max LL
    LLend = [paramEstim.LL];
    [LLMax idxMax] = max(LLend);

    %% draw
    clf
    for trial = 1:numParams
        subplot(331); hold on;

        semilogy(real(paramEstim(trial).LLs), 'color', lineColor);
        if HIGH_LIGHT_BEST_PARAM &(trial == idxMax)
            semilogy(real(paramEstim(trial).LLs), 'color', highColor);
        else
            semilogy(real(paramEstim(trial).LLs), 'color', lineColor);
        end
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

        plot(trial, paramEstim(trial).LLs(1), 'x', 'color', lineColor)
        if HIGH_LIGHT_BEST_PARAM &(trial == idxMax)
            plot(trial, paramEstim(trial).LL,'o', 'color', highColor);
        else
            plot(trial, paramEstim(trial).LL,'o', 'color', lineColor);
        end

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
    if ~isempty(paramInit)
        subplot(334); hold on
        plot([paramInit(:).alpha], 'xk', 'color', lineColor);
        subplot(335);  hold on
        plot([paramInit(:).beta], 'xk', 'color', lineColor);
        subplot(336); hold on
        plot([paramInit(:).Vrev], 'xk', 'color', lineColor);
    %     subplot(337); hold on
    %     plot(reshape([paramInit(:).Ke],[],numParams), 'xk--', 'color', lineColor);
        subplot(338); hold on
        plot([paramInit(:).Q], 'xk', 'color', lineColor);
        subplot(339); hold on
        plot([paramInit(:).R], 'xk', 'color', lineColor);
    end

    % estimated params
    subplot(334); hold on
    plot([paramEstim(:).alpha], 'o:', 'color', lineColor); title ('\alpha'); xlabel('trial');box off
    subplot(335); hold on
    plot([paramEstim(:).beta], 'o:', 'color', lineColor); title ('\beta'); xlabel('trial'); box off
    subplot(336); hold on
    plot([paramEstim(:).Vrev], 'o:', 'color', lineColor); title ('V_o'); xlabel('trial'); box off
    subplot(337); hold on
    idx = ~isnan([paramEstim(:).alpha]);
    KeHat = reshape([paramEstim(idx).Ke], nke,[]);
    %plot(reshape([paramEstim(:).Ke], nke,[]),'o-')
    if (size(KeHat,2)==1)
        plot(KeHat, 'o-', 'color', lineColor); 

    else
        errorbar(mean(KeHat'), std(KeHat'), 'color', lineColor); 

    end

    box off
    subplot(338); hold on
    plot([paramEstim(:).Q], 'o:', 'color', lineColor); title ('Q'); xlabel('trial'); box off
    subplot(339); hold on
    plot([paramEstim(:).R], 'o:', 'color', lineColor); title ('R'); xlabel('trial'); box off
    % subplot(259);
    % plot([paramEstim(:).Xo], 'o:'); title ('Xo'); xlabel('trial'); box off
    % subplot(2,5,10);
    % plot([paramEstim(:).Po], 'o:'); title ('Po'); xlabel('trial'); box off


    if (HIGH_LIGHT_BEST_PARAM)
        subplot(334); hold on
        plot(idxMax, paramEstim(idxMax).alpha, 'o:', 'color', highColor); title ('\alpha'); xlabel('trial');box off
        subplot(335); hold on
        plot(idxMax, paramEstim(idxMax).beta, 'o:', 'color', highColor); title ('\beta'); xlabel('trial'); box off
        subplot(336); hold on
        plot(idxMax, paramEstim(idxMax).Vrev, 'o:', 'color', highColor); title ('V_o'); xlabel('trial'); box off
        subplot(337); hold on
        plot(paramEstim(idxMax).Ke, '--', 'color', highColor);     box off
        subplot(338); hold on
        plot(idxMax, paramEstim(idxMax).Q, 'o:', 'color', highColor); title ('Q'); xlabel('trial'); box off
        subplot(339); hold on
        plot(idxMax, paramEstim(idxMax).R, 'o:', 'color', highColor); title ('R'); xlabel('trial'); box off
    end

    % plot true param
    if (~isempty(paramTrue))
        subplot(334);  hold on
        plot(repmat(paramTrue.alpha,numParams,1), '+r--'); 
        subplot(335); hold on
        plot(repmat(paramTrue.beta,numParams), '+r--'); 
        subplot(336); hold on
        plot(repmat(paramTrue.Vrev,numParams), '+r--'); 
        subplot(337); hold on
        plot(paramTrue.Ke, '+r--'); 
        subplot(338);  hold on
        plot(repmat(paramTrue.Q,numParams), '+r--'); 
        set(gca,'ylim', [0 max(paramTrue.Q*2, max([paramEstim(:).Q]))]);
        subplot(339); hold on
        plot(repmat(paramTrue.R,numParams), '+r--'); 
        set(gca,'ylim', [0 max(paramTrue.R*2, max([paramEstim(:).R]))]);
    % %     subplot(259); hold on
    % %     plot(repmat(paramTrue.Xo,numParams), '+r--'); 
    % %     subplot(2,5,10); hold on
    % %     plot(repmat(paramTrue.Po,numParams), '+r--'); 
    % %     set(gca,'ylim', [0 max(paramTrue.Po*2, max([paramEstim(:).Po]))]);
    end
else
    clf;
    [nn itr] = size(paramEstim);
    nke = nn-7;

    % log likelihood over EM-iterations
    subplot(331); hold on;
    semilogy(real(LLs));
    title (sprintf('L(\\Theta), final=%.3e',LLs(end)))


    subplot(332); hold on;
    semilogy(diff(real(LLs)));
    if length(LLs)>1
        title (sprintf('\\Delta L(\\Theta), final=%.3e',LLs(end)-LLs(end-1)))
    else
        title (sprintf('\\Delta L(\\Theta)'))
    end
    xlabel ('EM iteration')

    subplot(334);
    plot(paramEstim(1,1:itr), '.-');title ('\alpha');box off
    subplot(335);
    plot(paramEstim(2,1:itr), '.-');title ('\beta');box off
    subplot(336);
    plot(paramEstim(3,1:itr), '.-');title ('Vrev');box off
    subplot(337);
    plot(paramEstim(4:3+nke,1:itr)', '.-');title ('Ke');box off
    subplot(338);
    plot(paramEstim(4+nke,1:itr), '.-');title ('Q');box off
    xlabel ('EM iteration')
    subplot(339);
    plot(paramEstim(5+nke,1:itr), '.-');title ('R');box off
    xlabel ('EM iteration')

    % plot true param
    if (~isempty(paramTrue))
        subplot(334);
        hold on; plot(repmat(paramTrue.alpha,itr), 'r--');
        set(gca,'ylim', [min(0.9,min(paramEstim(1,1:itr))) 1])
        subplot(335);
        hold on; plot(repmat(paramTrue.beta,itr), 'r--');
        set(gca,'ylim', [0 paramTrue.beta*2])
        set(gca,'ylim', [min(paramTrue.beta*0.1,min(paramEstim(2,1:itr))) max(paramTrue.beta*1.1,max(paramEstim(2,1:itr)))])
        subplot(336);
        hold on; plot(repmat(paramTrue.Vrev,itr), 'r--');
        %set(gca,'ylim', [-0.01 0.01])
        set(gca,'ylim', [min(paramTrue.Vrev-0.01,min(paramEstim(3,1:itr))) max(paramTrue.Vrev+0.01,max(paramEstim(3,1:itr)))])
        subplot(337);
        hold on; plot(repmat(reshape(paramTrue.Ke,1,[]),itr,1), '--');
        subplot(338);
        hold on; plot(repmat(paramTrue.Q,itr), 'r--');
        set(gca,'ylim', [0 max(paramTrue.Q*2,max(paramEstim(4+nke,1:itr)))])
        xlabel ('EM iteration')
        subplot(339);
        hold on; plot(repmat(paramTrue.R,itr), 'r--');
        set(gca,'ylim', [0 max(paramTrue.R*2,max(paramEstim(5+nke,1:itr)))])
        xlabel ('EM iteration')
    % %     subplot(259);
    % %     hold on; plot(repmat(paramTrue.Xo,itr), 'r--');
    % %     set(gca,'ylim', [min(paramTrue.Xo-0.1,min(paramEstim(6+nke,1:itr))) max(paramTrue.Xo+0.1,max(paramEstim(6+nke,1:itr)))])
    % %     xlabel ('EM iteration')
    % %     subplot(2,5,10);
    % %     hold on; plot(repmat(paramTrue.Po,itr), 'r--');
    % %     set(gca,'ylim', [0 max(paramTrue.Po*2,max(paramEstim(7+nke,1:itr)))])
    % %     xlabel ('EM iteration')    


    end
end



% set(gcf, 'paperposition', [0 0 15 6]) 
% set(gcf, 'papersize', [15 6]) 
% saveas(1, sprintf('params_EM.pdf'));