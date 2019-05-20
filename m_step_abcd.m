function [A B C D Q R LL] = m_step_abcd(Y, U, Xs, Ps, Pcs, A, B, C, D, Q, R)


N = length(Y);

[dimX,~] = size(A);

    % for updating A and B
    S0 = sum(Xs(:,1:end-1),2);    % sum mean(Xn) from 1 to end-1
    %S1 = S0+Xs(:,end)-Xs(:,1);  % sum mean(Xn) from 2 to end
    

        
    % calc square sums 
    SS00 = 0;                   
    SS10 = 0;
    for i=1:N-1
        SS00 = SS00 + Xs(:,i)*Xs(:,i)' + Ps(:,:,i); %sum(Xn*Xn') from 1 to end-1
        SS10 = SS10 + Xs(:,i+1)*Xs(:,i)' +  Pcs(:,:,i+1); %sum(X{n+1}*Xn') from 1 to end-1
    end
    %sum(Xn*Xn') from 1 to end
    SS = SS00 + (Xs(:,end)*Xs(:,end)'+ Ps(:,:,end));
    %sum(Xn*Xn') from 2 to end
    SS11 = SS - (Xs(:,1)*Xs(:,1)'- Ps(:,:,1));   

    
    % calc terms that involves input
    % calc square sums 
    SXU00 = 0;                   
    SXU10 = 0;                   
    SUU00 = 0;
    for i=1:N-1
        SXU00 = SXU00 + Xs(:,i)*U(:,i)'; %sum(Xn*Un') from 1 to end-1
        SXU10 = SXU10 + Xs(:,i+1)*U(:,i)'; %sum(X{n+1}*Un') from 1 to end-1
        SUU00 = SUU00 + U(:,i)*U(:,i)';   %sum(Un*Un') from 1 to end-1
    end
    SXU = SXU00 + Xs(:,end)*U(:,end)';
    SUU = SUU00 + U(:,end)*U(:,end)';


    % calc terms that involves measurement (for updating D and R)
    SYU = 0;
    SYY = 0;
    SYX = 0;
    for i=1:N
        SYU = SYU + Y(:,i)*U(:,i)';
        SYY = SYY + Y(:,i)*Y(:,i)';
        SYX = SYX + Y(:,i)*Xs(:,i)';
    end

    
    % calc quadratic terms
    quad1 = SS11+A*SS00*A'+B*SUU00*B'  - 2*SS10*A' + 2*A*SXU00*B' -2*SXU10*B' ;
    quad2 = SYY + C*SS*C' + D*SUU*D'- D*SYU' - SYU*D' - C*SYX' -SYX*C'  ...
             + C*SXU*D' + D*SXU'*C';
         

    % update A and B
    AB = [SS10 SXU10]/[SS00 SXU00; SXU00' SUU00];
    if ~sum(isnan(AB(:)))
%         if (updateParam(1))
            A = AB(:,1:dimX);
%         end
%         if (updateParam(2))
            B = AB(:,dimX+1:end);
%         end
    end

    % update D
%     if (updateParam(4))
%         D = (SYU-C*SXU)/SUU;
%     end

    % update Q    
%     if (updateParam(5))
        Q = quad1 / (N-1);
        Q = diag(diag(Q));       % diagonalize

%     end

    % update R
%     if (updateParam(6))
        R = (quad2)/N;
        R = diag(diag(R));        % diagonalize
%     end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% calc LL & check convergence 
    L1 =  trace(Q\quad1) + (N-1)*log(abs(det(Q)));
    L2 =  trace(R\quad2) + N*log(abs(det(R)));
    
    LL = -0.5*(L1 + L2  );