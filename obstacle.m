function [uh,in,ierr] = obstacle(psi,g,f,tol,fd,h0,p,t,varargin);
%OBSTACLE  Solve the obstacle problem
%   u >= psi in D, - Laplacian u = f on {u>psi}, u=g on boundary of D,  
%using triangulation described by  p,t.  Assumes psi>= g on boundary of D.
%Reports solution in  uh, indices of interior points in  in, and error at
%each iteration in  ierr.
%Example:
%  >> fd=@(p) sqrt(sum(p.^2,2))-1;  psi=@(p) -sum((3*p).^4,2)+1;  f=@(p) 0;
%  >> h0=0.06; figure(1), [p,t]=distmesh2d(fd,@huniform,h0,[-1,-1;1,1],[]);
%  >> [uh,in,ierr]=obstacle(psi,f,f,1e-6,fd,h0,p,t);  
%  >> figure(2), semilogy(1:length(ierr),ierr,'o'), xlabel j, ylabel('iter err');
%  (old style:
%     fd=inline('sqrt(sum(p.^2,2))-1','p');  
%     psi=inline('-sum((3*p).^4,2)+1','p');  f=inline('0','p');
%     h0=0.06; figure(1), [p,t]=distmesh2d(fd,@huniform,h0,[-1,-1;1,1],[]);
%     [uh,in,ierr]=obstacle(psi,f,f,1e-6,fd,h0,p,t);  
%     figure(2), semilogy(1:length(ierr),ierr,'o'), xlabel j, ylabel('iter err');
%
%   See also: POISSONDN.
%ELB 12/3/04

% square peg:  psi=@(p) 2*((p(:,1).^2<.1)&(p(:,2).^2<.1))-1;

omega=1.75;     % found by trial and error
maxiter=300;

% use poissonDN to get unconstrained stiffness, load
[uh,in,A,b]=poissonDN(f,g,@(p)(0),fd,@(p)(-1),h0,p,t,varargin{:});
U=triu(A,1); L=tril(A,-1); d=diag(A);     % U, L sparse
if any(d==0), error('stiffness matrix has zero on diagonal'), end;

% first guess is max(uh,psi)
N=sum(in>0);  ps=zeros(N,1);
for j=1:N,  ps(j)=feval(psi, p(find(in==j),:) ); end   
uold=max(uh(in>0),ps);  unew=uold;  omcomp=1-omega;  ierr=[];

% iterate: constrained point over-relaxation
for l=1:maxiter+1
    Ux=U*uold;
    for j=1:N
        utemp=(b(j)-L(j,1:j-1)*unew(1:j-1)-Ux(j))/d(j);  % Gauss-Seidel
        unew(j)=max(omcomp*uold(j)+omega*utemp,ps(j));  end
    er=max(abs(unew-uold));  ierr=[ierr er];
    if er<tol, break, end
    if l>maxiter, warning('max number of iterations reached'), break, end 
    uold=unew;  end

uh(in>0)=unew;
h=trimesh(t,p(:,1),p(:,2),uh);  set(h,'FaceAlpha',0.3)  % plot transparent
xy=[get(gca,'Xlim') get(gca,'YLim')];  hold on;
trisurf(t,p(:,1),p(:,2),psi(p),uh);  axis([xy min(uh) max(uh)]);  hold off;
