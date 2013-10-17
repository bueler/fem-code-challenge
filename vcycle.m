function [uu] = vcycle(lev, A, u, b, R, P, nu1, nu2)
%VCYCLE Apply one V-cycle from level=lev using nu1 presmoothing and
%       nu2 postsmoothing steps. A,u,b,R,P are cell arrays containing
%	matrices on each level. See fmg.m for usage examples.
%
%   See also: FMG.
%Jed 2004-11-16

uu = smooth(A{lev}, u{lev}, b{lev}, nu1);
res = b{lev} - A{lev}*uu;
b{lev-1} = R{lev-1}*res;
if (lev == 2)
	delta = A{lev-1} \ b{lev-1};
else
	u{lev-1} = zeros(size(u{lev-1}));
	delta = vcycle(lev-1, A, u, b, R, P, nu1, nu2);
end
uu = uu + P{lev-1}*delta;
uu = smooth(A{lev}, uu, b{lev}, nu2);


function uu = smooth(A, u, b, nu)
%SMOOTH		Apply nu iterations of Gauss-Seidel to the system Au=b
%Jed 2004-11-16;  ELB 11/17/04

L = tril(A);  U = triu(A,1);  uu = u;
for i=1:nu
	uu = L \ (b - U*uu);
end
