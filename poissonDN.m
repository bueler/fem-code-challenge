function [uh,un,A,b]=poissonDN(f,gD,gN,fd,fGam,h0,p,t,varargin)
%POISSONDN  Solve Poisson's equation on (open) domain D by the FE method:
%    -Laplacian u = f on D
%with nonhomogeneous Dirichlet and Neumann boundary conditions
%    u=gD  on  GamD    and    du/dn=gN  on  GamN
%where GamD is a closed subset of the boundary and GamN=bdry(D)-GamD.  Uses 
%triangulation of D  given by p,t, a mesh size h0, and a signed distance 
%function fd.  Note fd(x)<0 iff x is in D.  Similarly, fGam(x)<=0 if x in 
%GamD and fGam(x)>0 if x in GamN.  Requires gD(x) defined for all x in 
%GamD and gN(x) defined for all x in GamN.  Note f(x) should be defined for
%all nodes.  Returns an approximate solution uh defined on all nodes and 
%the numbered unknown nodes in un.  Unknown nodes are those in interior 
%and those in GamN.
%
%Example: To solve  -lap u = 0  on  D=[0,1]x[0,1]  and u=x^2-y^2 on bdry D,
%  >> fd=@(p) drectangle(p,0,1,0,1);  gD=@(p)(p(:,1).^2-p(:,2).^2);
%  >> [p,t]=distmesh2d(fd,@huniform,.1,[0,0;1,1],[0,0;0,1;1,0;1,1]);
%  >> uh=poissonDN(@(p)(0),gD,@(p)(0),fd,@(p)(-1),.1,p,t);
%  >> err=max(abs(uh-gD(p)))
%  (old style: 
%     fd=inline('drectangle(p,0,1,0,1)','p');  fGam=inline('-1','p');
%     f=inline('0','p');  gD=inline('p(:,1).^2-p(:,2).^2','p');
%     [p,t]=distmesh2d(fd,@huniform,.1,[0,0;1,1],[0,0;0,1;1,0;1,1]);
%     uh=poissonDN(f,gD,f,fd,fGam,.1,p,t);
%     err=max(abs(uh-gD(p)))
%Further convergence examples in  TESTDNSQ, TESTDNTENT, TESTDNDISC.
%
%   See also: POISSONV2, DISTMESH2D, NUMBERMESH, TRIMESH.
%ELB 11/19/04

geps=.001*h0;  int=(feval(fd,p,varargin{:}) < -geps);    % int true if node in interior
inGamD=(feval(fGam,p) < +geps)&(~int);       % inGamD true if node in GamD
inGamN=~(int|inGamD);                        % inGamN true if node in GamN
if sum(inGamD)==0, error('no unique soln; needs Dirichlet points'), end
Np=size(p,1);  uh=zeros(Np,1);  ff=uh;  un=uh;  % Np=total # of nodes
N=sum(~inGamD);  un(~inGamD)=(1:N)';            % N=# of unknowns
for j=1:Np % eval f once for each node; fill in known bdry vals on Gam1
    ff(j)=feval(f,p(j,:)); 
    if inGamD(j), uh(j)=feval(gD,p(j,:)); end, end   

% loop over triangles to set up stiffness matrix A and load vector b
A=sparse(N,N);  b=zeros(N,1);
for n=1:size(t,1)
    j=t(n,1);  k=t(n,2);  l=t(n,3);  vj=un(j);  vk=un(k);  vl=un(l);
    J=[p(k,1)-p(j,1), p(l,1)-p(j,1); p(k,2)-p(j,2), p(l,2)-p(j,2)];
    ar=abs(det(J))/2;  C=ar/12;  Q=inv(J'*J);  fT=[ff(j) ff(k) ff(l)];
    % go through nodes and compute stiffness and Dirichlet contribution
    if vj>0
        A(vj,vj)=A(vj,vj)+ar*sum(sum(Q)); b(vj)=b(vj)+C*fT*[2 1 1]';
        if vk>0
            A(vj,vk)=A(vj,vk)-ar*sum(Q(:,1));  A(vk,vj)=A(vj,vk); end
        if vl>0
            A(vj,vl)=A(vj,vl)-ar*sum(Q(:,2));  A(vl,vj)=A(vj,vl); end
    else % pj in GamD
        if vk>0, b(vk)=b(vk)+uh(j)*ar*sum(Q(:,1)); end
        if vl>0, b(vl)=b(vl)+uh(j)*ar*sum(Q(:,2)); end, end
    if vk>0
        A(vk,vk)=A(vk,vk)+ar*Q(1,1);  b(vk)=b(vk)+C*fT*[1 2 1]';
        if vl>0
            A(vk,vl)=A(vk,vl)+ar*Q(1,2);  A(vl,vk)=A(vk,vl); end
    else % pk in GamD
        if vj>0, b(vj)=b(vj)+uh(k)*ar*sum(Q(:,1)); end
        if vl>0, b(vl)=b(vl)-uh(k)*ar*Q(1,2); end, end
    if vl>0
        A(vl,vl)=A(vl,vl)+ar*Q(2,2);  b(vl)=b(vl)+C*fT*[1 1 2]';
    else % pl in Gam1
        if vj>0, b(vj)=b(vj)+uh(l)*ar*sum(Q(:,2)); end
        if vk>0, b(vk)=b(vk)-uh(l)*ar*Q(1,2); end, end
    % now add Neumann contribution
    if inGamN(j)&&inGamN(k)&&inGamN(l)
        error('illegal triangle with all vertices in Neumann bdry'), end
    if inGamN(j)
        if ~int(k), b(vj)=b(vj)+Ncntrb(gN,p(j,:),p(k,:)); end
        if ~int(l), b(vj)=b(vj)+Ncntrb(gN,p(j,:),p(l,:)); end,  end
    if inGamN(k)
        if ~int(j), b(vk)=b(vk)+Ncntrb(gN,p(k,:),p(j,:)); end
        if ~int(l), b(vk)=b(vk)+Ncntrb(gN,p(k,:),p(l,:)); end,  end
    if inGamN(l)
        if ~int(j), b(vl)=b(vl)+Ncntrb(gN,p(l,:),p(j,:)); end
        if ~int(k), b(vl)=b(vl)+Ncntrb(gN,p(l,:),p(k,:)); end,  end
end

uh(~inGamD)=A\b;                             % solve for FE solution
trimesh(t,p(:,1),p(:,2),uh), axis tight      % display

function w=Ncntrb(gN,p,q); % compute Neumann contribution by Simpson's rule
w=norm(p-q)*( feval(gN,p) + 2*feval(gN,(p+q)/2) )/6;
