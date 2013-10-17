function err=testDNsq(h0)
%TESTDNSQ  Test poissonDN.m on a square with purely Dirichlet,
%nonhomogeneous boundary values.  In particular, 
%  D=[0,1] x [0,1];   - lap u = 0;   u=x^2-y^2 on bdry D
%which has solution  u=x^2-y^2.  Returns maximum node error.
%
%Example:
%  >> m=6;  h0=.5*(3/5).^(0:m),  err=testDNsq(h0)
%has execution time of a couple of minutes.
%
%  See also:  POISSONDN, DISTMESH2D, TESTDNSQMIX, TESTDNTENT
%ELB 11/8/04

fd=inline('drectangle(p,0,1,0,1)','p');  fGam=inline('-1','p');
f=inline('0','p');  gD=inline('p(:,1).^2-p(:,2).^2','p');

clear err,  figure(1),  M=length(h0);
for j=1:M
    disp(['h0       = ' num2str(h0(j))]),  clf,  tic
    [p,t]=distmesh2d(fd,@huniform,h0(j),[0,0;1,1],[0,0;0,1;1,0;1,1]);
    clf,  [uh,un]=poissonDN(f,gD,f,fd,fGam,h0(j),p,t);  drawnow
    err(j)=max(abs(uh-gD(p)));
    disp(['time     = ' num2str(toc)])
    disp(['node err = ' num2str(err(j))])
end

figure(2), loglog(h0,err,'o',h0,err(M)*(h0/h0(M)).^2,'--')
grid on, axis tight,  title('max node err;  dashed line is  C h_0^2')
