%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Function: bdyrefine.m
%%
%% Purpose:  Refines a triangulation using edge bisection while also
%%           adjusting boundary points to lie on a domain boundary.
%%
%% Input:  p:     matrix of points (nx2)
%%         t:     matrix of indices into p defining a triangulation (Nx3)
%%         fd:    distance function a-la distmesh2d
%%
%%
%% Output: rp:    an augmented list of points in a refined triangluation.
%%                The first N points will be the same as those in p,
%%                possibly shifted at the boundary
%%         rt:    a list of triangles in the refined triangulation
%%         e:     Edgelist used for refinement
%%         ind:   Logical mask for interior elements of p
%%
%% Example:
%%
%%  fd=inline('sqrt(sum(p.^2,2))-1','p')
%%  [p,t] = distmesh2d( fd, @huniform, 0.15, [-1,-1;1,1],[]);
%%  [rp,rt] = bdyrefine( p, t, fd );
%%
%% if edgelist and logical mask for interior points are desired:
%%  [rp,rt,e,ind] = bdyrefine( p, t, fd );
%%
%% DAM 10/26
%% Modified: Jed 2004-11-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rp,rt,e,ind] = bdyrefine(p,t,fd,varargin)

  [e,te] = edgelist(p,t);

  N = size( p, 1);
  rp = zeros( size(p,1)+size(e,1), 2 );
  rp( 1:size(p,1),:) = p;

  % Add new boundary points to the end of the point list
  rp(N+1:end,:) = 0.5*(p(e(:,1),:) + p(e(:,2),:));

  % Find boundary edges by traversing edges and counting
  bdy_count = zeros( size(e,1), 1 );
  for( i=1:size(te,1) )
    for( j=1:3 )
      bdy_count( te(i,j) ) = bdy_count( te(i,j) ) + 1;
    end
  end

  % A boundary point on the fine grid is either an endpoint of a coarse
  % grid boundary edge or a midpoint of a coarse boundary edge.
  bdy_point=logical(zeros( size( rp, 1 ), 1 ));
  bdy_edge = bdy_count == 1;
  bdy_point( e(bdy_edge,1) ) = (1==1);
  bdy_point( e(bdy_edge,2) ) = (1==1);
  bdy_point( N+1:size(bdy_point,1)) = bdy_edge;

  % Construct 4 triangles for every triangle in the original mesh
  rt = zeros( 4*size(t,1), 3 );
  for i=1:size(t,1)
    rt(4*(i-1)+1,:) = [ N+te(i,3), t(i,1), N+te(i,1) ];
    rt(4*(i-1)+2,:) = [ N+te(i,1), t(i,2), N+te(i,2) ];
    rt(4*(i-1)+3,:) = [ N+te(i,2), t(i,3), N+te(i,3) ];
	rt(4*(i-1)+4,:) = [ N+te(i,1), N+te(i,2), N+te(i,3) ];
  end

  % Use numerical gradient technique stolen from distmesh2d to project
  % boundary points onto the boundary.
  deps=sqrt(eps)/sqrt(size(rt,1))*.001;
  d=feval(fd,rp(bdy_point,:),varargin{:});
  dgradx=(feval(fd,[rp(bdy_point,1)+deps,rp(bdy_point,2)],varargin{:})-d)/deps;
  dgrady=(feval(fd,[rp(bdy_point,1),rp(bdy_point,2)+deps],varargin{:})-d)/deps;
  rp(bdy_point,:)=rp(bdy_point,:)-[d.*dgradx,d.*dgrady];
  ind = (bdy_count == 2);
