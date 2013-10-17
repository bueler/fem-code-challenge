function err=testDNtent(h0)
%TESTDNTENT  Test of poissonDN.m on a rectangle with nonhomogeneous 
%Dirichlet boundary values but where the exact solution is not in
%H^2(D).  In particular we solve Laplace's equation  - lap u = 0  on
%D=[0,pi] x [0,1]  and
%                                            / x,      0 <= x <= pi/2,
%   u(x,0) = u(0,y) = u(pi,y) = 0,  u(x,1)= |
%                                            \ pi - x, pi/2 <= x <= pi.
%Compares to the solution found by separation of variables.  Returns max
%node err.
%
%Example: (execution time several minutes)
%  >> m=5; h0=.5*(3/5).^(0:m); err=testDNtent(h0)
%
%  See also:  POISSONDN, DISTMESH2D, TESTDNSQ, TESTDNSQMIX
%ELB 11/8/04

global h
fd=inline('drectangle(p,0,pi,0,1)','p');  fGam=inline('-1','p');
f=inline('0','p');
figure(1), clf, M=length(h0);
for j=1:M
    h=h0(j);  disp(['h0       = ' num2str(h)]), tic, clf
    [p,t]=distmesh2d(fd,@huniform,h,[0,0;pi,1],[0,0;0,1;pi,0;pi,1]);
    clf, [uh,un]=poissonDN(f,@gD,f,fd,fGam,h,p,t);  drawnow
    err(j)=max(abs(uh-getexact(500,p)));
    disp(['time     = ' num2str(toc)])
    disp(['node err = ' num2str(err(j))])
end
figure(2), clf, loglog(h0,err,'o',h0,err(M)*(h0/h0(M)).^2,'--')
grid on, title('node err;  dashed line is C (h_0)^2')

function z=gD(p);
global h
z=zeros(size(p,1),1);  top=( abs(p(:,2)-1)<.001*h );
lt= top && (p(:,1) < pi/2);   rt= top && (p(:,1) >= pi/2);
z(lt)=p(lt,1); z(rt)=pi-p(rt,1);
    
function u=getexact(N,p);
k=1:N; c=4*sin(k*pi/2)./(pi*k.^2);  u=zeros(size(p,1),1); 
for k=N:-1:1
    u=u+c(k)*sin(k*p(:,1)).*(sinh(k*p(:,2))/sinh(k));
end
