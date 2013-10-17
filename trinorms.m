function [L2,semi1,H1]=trinorms(v,p,t);
%TRINORMS  A utility to compute the L^2 norm, the L^2 norm of the gradient 
%(a seminorm) and the H^1 norm of a function on a triangular mesh.  Note 
%these quantities are related by H1=sqrt(L2^2+semi1^2).
%
%Example:  To compute the norms of piecewise linear approximations to 
%v(x,y)=xy if D is the square centered at the origin and with side 2,
% >> fd=inline('drectangle(p,-1,1,-1,1)','p');
% >> [p,t]=distmesh2d(fd,@huniform,0.1,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
% >> [L2,semi1,H1]=trinorms(p(:,1).*p(:,2),p,t)
%The exact L^2 norm of v is 2/3.  The exact value of the L^2 norm of the 
%gradient of v is sqrt(8/3)=1.63299316185545.  The exact value of the H^1 
%norm of v is 2*sqrt(7)/3=1.76383420737639.
%
%   See also: DISTMESH2D, TRIMESH, POISSONV2.
%ELB 10/31/04

L2=0;  semi1=0;  M=[2 1 1; 1 2 1; 1 1 2];
for n=1:size(t,1)
    j=t(n,1); k=t(n,2); l=t(n,3);  y=[v(j) v(k) v(l)];
    J=[p(k,1)-p(j,1), p(l,1)-p(j,1); p(k,2)-p(j,2), p(l,2)-p(j,2)];
    Jdet=abs(det(J));
    L2=L2+(Jdet/24)*(y*(M*y'));
    Q=inv(J'*J);  dy=v(j)*[-1 -1]+v(k)*[1 0]+v(l)*[0 1];
    semi1=semi1+(Jdet/2)*dy*Q*dy';
end

% at this point L2, semi1 are squares of (semi)norms
H1=sqrt(semi1+L2);  L2=sqrt(L2);  semi1=sqrt(semi1);
