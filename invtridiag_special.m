function [L0,L1] = invtridiag_special(M11,M21,M22,nd)
% [L0,L1] = invtridiag_special(M11,M21,M22,nd)
%
% Find the central ndiags of the inverse of a tri-diagonal matrix using
% method of Rybicki & Hummer (1991).  Uses special non-recursive formula
% developed by Memming.
%
% Special format: assumes we have tridiagonal of the form:
%
% [A B 0 0 0
%  B C B 0 0
%  0 B C B 0
%     . . .
%  0 0 B C B
%  0 0 0 B A]
%
% Inputs: 
%   M11 - M(1,1) = top left entry
%   M21 - M(2,1) = 1st off-diagonal term
%   M22 - M(2,2) = 2nd off-diagonal term
%   nd - number of dimensions of full(nd x nd) matrix
%
% Output: L0 - main diagonal
%         L0 - off-diagonal

a = M11/M21;
q = M22/M21;

sq = sqrt(q^2-4);
rp = (q+sq)./(q-sq);
rr = rp.^(1:(nd-2));

d = ((a*(sq-q)+2)+(a*(sq+q)-2)*rr)./ ...
    (((sq+q)-2*a)+((sq-q)+2*a)*rr);
d = [inf a d inf]';

% Now form L0 itself
de = 1./(d.*flipud(d));
L0 = 1./((1-de(2:nd+1)).*(M22-M21./d(1:nd)));
L0(1) = L0(1)/(M11-M21/d(1))*(M22-M21/d(1));
L0(nd) = L0(1);

% Form off-diagonal term
L1 = -L0(2:nd)./d(2:nd);

