function err=testDNsqmix(h0)
%TESTDNSQMIX  Test poissonDN.m on a square with mixed Dirichlet and 
%Neumann nonhomogeneous boundary values.  In particular, 
%  D=[0,1] x [0,1];   - lap u = 0;   
%  du/dn = - 2  on  GamN = { (x,1) | 0 < x < 1 };
%  u=x^2-y^2 on GamD = bdry(D) - GamN
%which has solution  u=x^2-y^2.
%
%Example:
%  >> m=6;  h0=.5*(3/5).^(0:m),  err=testDNsqmix(h0)
%has execution time of a couple of minutes.
%
%  See also:  POISSONDN, DISTMESH2D, TESTDNSQ, TESTDNTENT
%ELB 11/8/04

fd=inline('drectangle(p,0,1,0,1)','p');
fGam=inline('.25 - (p(:,1)-.5).^2 - (p(:,2)-1).^2','p');
f=inline('0','p');  
gD=inline('p(:,1).^2-p(:,2).^2','p');  gN=inline('-2*p(:,2)','p');

clear err,  M=length(h0);
for j=1:M
    disp(['h0       = ' num2str(h0(j))]),  figure(1),  clf,  tic
    [p,t]=distmesh2d(fd,@huniform,h0(j),[0,0;1,1],[0,0;0,1;1,0;1,1]);
    figure(2),  clf
    [uh,un]=poissonDN(f,gD,gN,fd,fGam,h0(j),p,t);  drawnow
    err(j)=max(abs(uh-gD(p)));
    disp(['time     = ' num2str(toc)])
    disp(['node err = ' num2str(err(j))])
end

figure(3), loglog(h0,err,'o',h0,err(M)*h0/h0(M),'--')
grid on, axis tight,  title('max node err;  dashed line is  C (h_0)^1')
