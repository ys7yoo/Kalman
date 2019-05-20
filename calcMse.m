function  [mse dcDiff] = calcMse(X, Xtrue, subtractDC)
if nargin <3
    subtractDC=0;
end
   
if subtractDC
    Xdc=mean(X);
    XtrueDc=mean(Xtrue);
    
    X = bsxfun(@minus,X,Xdc);
    Xtrue = bsxfun(@minus,Xtrue,XtrueDc);

    
    dcDiff = Xdc-XtrueDc;
end
    
% calc mse of Ke 
M = size(Xtrue,1);
Err = bsxfun(@minus,X,Xtrue);
mse = diag(Err'*Err)' / M;

