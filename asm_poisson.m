function [A, b] = asm_poisson(f,p,t,ind);
%ASM_POISSON  Assemble for Poisson's equation on a domain D by the FE method:
%   - Laplacian u = f on D, u=0 on boundary of D
%using triangulation described by  p,  t,  and indices of unknown points 
%ind.  Returns stiffness matrix  A  and load vector  b.  Modification of 
%poissonv2  useful in multigrid.
%Example:  (see  domg  and  fmg  for nontrivial use):
%  >> fd=@(p) sqrt(sum(p.^2,2))-1;
%  >> [p,t] = distmesh2d(fd,@huniform,0.5,[-1,-1;1,1],[]);
%  >> [rp,rt,e,ind] = bdyrefine(p,t,fd);
%  >> [A,b]=asm_poisson(@(p)(4),rp,rt,ind);  spy(A)
%  (old style:
%     fd=inline('sqrt(sum(p.^2,2))-1','p')
%     [p,t] = distmesh2d(fd,@huniform,0.5,[-1,-1;1,1],[]);
%     [rp,rt,e,ind] = bdyrefine(p,t,fd);
%     [A,b]=asm_poisson(inline('4','p'),rp,rt,ind);  spy(A)
%
%   See also: POISSONV2, FMG, DOMG, BDYREFINE.
%ELB 10/15/04
%Modified: Jed 2004-11-16
%Comments modified: ELB 2004-11-17

Np=size(p,1);  N=sum(ind);        % Np=# of nodes;  N=# of interior nodes
in=zeros(Np,1);  in(ind)=(1:N)';  % number the interior nodes
for j=1:Np, ff(j)=feval(f,p(j,1:2)); end   % eval f once for each node

% loop over triangles to set up stiffness matrix A and load vector b
A=sparse(N,N);  b=zeros(N,1);
for n=1:size(t,1)
    j=t(n,1); k=t(n,2); l=t(n,3); J=in(j); K=in(k); L=in(l);
    Jac=[p(k,1)-p(j,1), p(l,1)-p(j,1); p(k,2)-p(j,2), p(l,2)-p(j,2)];
    ar=abs(det(Jac))/2;  C=ar/12;  Q=inv(Jac'*Jac);  fT=[ff(j) ff(k) ff(l)];
    if J>0
        A(J,J)=A(J,J)+ar*sum(sum(Q));  b(J)=b(J)+C*fT*[2 1 1]';
        if K>0, A(J,K)=A(J,K)-ar*sum(Q(:,1));  A(K,J)=A(J,K); end
        if L>0, A(J,L)=A(J,L)-ar*sum(Q(:,2));  A(L,J)=A(J,L); end, end
    if K>0
        A(K,K)=A(K,K)+ar*Q(1,1);  b(K)=b(K)+C*fT*[1 2 1]';
        if L>0, A(K,L)=A(K,L)+ar*Q(1,2);  A(L,K)=A(K,L); end, end
    if L>0, A(L,L)=A(L,L)+ar*Q(2,2);  b(L)=b(L)+C*fT*[1 1 2]'; end
end
