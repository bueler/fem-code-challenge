%DOMG	Use Full Multigrid (FMG) to solve Poisson's equation on a disc in 
%    case where exact solution is known.  Record times and errors of FMG 
%    versus direct ("A\b") solution of sparse FEM equations.  (Only 
%    intended as an example of use of fmg.m and asm_poisson.m.)
%
%   See also: FMG, ASM_POISSON.
%Jed 2004-11-16


nlev = 7; % default to 6; if  nlev==6  then a couple of minutes
h0 = 0.4; 

f=inline('4','p');  fd=inline('sqrt(sum(p.^2,2))-1','p');
[p,t] = distmesh2d_noplot(fd, @huniform, h0, [-1,-1;1,1],[]);
time = zeros(nlev,5); err = zeros(nlev,2);
for lev=2:nlev
	lev
	[uh,pp,tt,A,b,ind,time(lev,1:4)] = fmg(@asm_poisson,f,fd,h0,p,t,lev,1,1,4);
	u=1-sum(pp.^2,2); % Exact solution
	err(lev,1) = max(abs(uh-u));	% Error of FMG solution
	tic, up = zeros(size(pp,1),1); up(ind) = A\b; % FEM direct solution
	time(lev,5) = toc; err(lev,2) = max(abs(up-u)); % Error of direct solution
end
disp('format of "time": Refinement, assembly, MG prep, FMG solve, A\b solve')
time=time(2:end,:), 

disp('format of "err": FMG err, A\b err')
err=err(2:end,:)
