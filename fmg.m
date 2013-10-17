function [u_fine, p_fine, t_fine, A, b, ind, time] = ...
    fmg(assemble, f, fd, h0, p, t, nlev, nu1, nu2, mu);
%FMG	Full Multigrid Method
% Inputs:
%	assemble	Function handle for   [A b]=assemble(f,p,t,ind) 
%             assemble creates stiffness matrix A and load vector b for
%             PDE problem w. non-homogeneous term f on a triangulation p,t
%             and interior points p(ind).
%	f		Non-homogeneous term.
%	fd		Signed distance function (see Persson & Strang's distmesh2d.m).
%	h0		Coarse mesh size
%	p		Point list for coarse triangulation
%	t		Triangle list for coarse triangulation
%	nlev	Number of levels of refinement to use (>=2)
%	nu1		Number of presmoothing iterations for V-cycle
%	nu2		Number of postsmoothing iterations for V-cycle
%	mu		Number of V-cycles
% Outputs:
%	u_fine  Value of approximate solution on finest grid
%	p_fine	Point list for finest triangulation
%	t_fine	Triangle list for finest triangulation
%	A		Stiffness matrix on finest triangulation
%	b		Load vector on finest triangulation
%	ind		Logical mask for interior points of p: u(ind)=A\b
%	time	1x4 Time vector: Refinement, assembly, MG prep, FMG solve
%
%   See also: ASM_POISSON, VCYCLE, DOMG, DISTMESH2D.
%Jed 2004-11-16

tic; time=zeros(1,4);
AA = cell(nlev, 1); uu = AA; bb = AA; hh = AA; NN = AA;
RR = cell(nlev-1, 1); PP = RR; u = AA;
hh{1} = h0; ind = feval(fd,p) < -.001*h0; n = sum(ind);
for lev = 2:nlev
	[rp,rt,e,int] = bdyrefine(p,t,fd);
	ind = [ind;int];
	N = size(p,1); NN{lev-1} = N; hh{lev} = hh{lev-1}/2;
	nc = n; n = sum(ind);
	in = zeros(size(ind)); in(ind) = 1:n;
	j = reshape(repmat(N+1:N+size(e,1),[2 1]),[2*size(e,1) 1]);
	kj = [reshape(e',[prod(size(e)) 1]) j];
	kj = in(kj); kj = kj(kj(:,1) ~= 0, :);
	P = speye(n,nc); P(sub2ind(size(P),kj(:,2),kj(:,1))) = 0.5;
	PP{lev-1} = P; RR{lev-1} = speye(nc,n);
	p = rp; t = rt;
end
NN{nlev} = size(p,1);
time(1)=toc; tic	% Done Refining
[A, b] = feval(assemble,f,rp,rt,ind);
time(2)=toc; tic	% Assembly complete
AA{nlev} = A; bb{nlev} = b;
for lev = nlev-1:-1:1
	AA{lev} = RR{lev} * (AA{lev+1} * PP{lev});
	bb{lev} = 2 * RR{lev} * bb{lev+1};
end
time(3)=toc; tic	% MG Prep complete
for lev = 1:nlev, u{lev} = zeros(NN{lev},1); end
uu{1} = AA{1} \ bb{1};
for lev = 2:nlev
	uu{lev} = PP{lev-1} * uu{lev-1};
	for s=1:mu
		uu{lev} = vcycle(lev, AA, uu, bb, RR, PP, nu1, nu2);
	end
end
time(4) = toc;		% FMG Solve complete
u{nlev}(ind) = uu{nlev};
%trisurf(rt,rp(:,1),rp(:,2),u{nlev})
u_fine = u{nlev}; p_fine = p; A=AA{nlev}; b = bb{nlev}; t_fine = t;
