%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%    File: refine.m
%%
%% Purpose: Refines a triangulation using edge bisection.
%%
%% Input:  p:     matrix of points (nx2)
%%         t:     matrix of indices into p defining a triangulation (Nx3)
%%
%% Output: rp:    an augmented list of points in a refined triangluation.  
%%                The first N points will be the same as those in p.
%%         rt:    a list of triangles in the refined triangulation
%%
%%
%% Example: 
%%
%%  fd=inline('sqrt(sum(p.^2,2))-1','p');
%%  [p,t] = distmesh2d( fd, @huniform, 0.15, [-1,-1;1,1],[]);
%%  [rp,rt] = refine( p, t ); 
%%
%% DAM 10/26
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rp, rt] = refine( p, t )

  [e, te ] = edgelist( p, t );
  N = size( p, 1);
  rp = zeros( size(p,1)+size(e,1), 2 );
  rp( 1:size(p,1),:) = p;
  for i=1:size(e,1)
    rp(N+i,:) = (p(e(i,1),:)+p(e(i,2),:))/2;
  end
  
  rt = zeros( 4*size(t,1), 3 );  
  for i=1:size(t,1)
    rt(4*(i-1)+1,:) = [ N+te(i,3), t(i,1), N+te(i,1) ];
    rt(4*(i-1)+2,:) = [ N+te(i,1), t(i,2), N+te(i,2) ];
    rt(4*(i-1)+3,:) = [ N+te(i,2), t(i,3), N+te(i,3) ];
	rt(4*(i-1)+4,:) = [ N+te(i,1), N+te(i,2), N+te(i,3) ];
  end
