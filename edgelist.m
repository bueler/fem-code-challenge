%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%    File: edgelist.m
%%
%% Purpose: Compute a list of edges in a triangulation
%%
%% Input:  p:     list of points (nx2)
%%         t:     list indices into p defining a triangulation (Nx3)
%%
%% Output: e:     for each edge a list of indices into p defining the
%%                endpoints of the edge (mx2)
%%         te:    for each triangle in the triangluation a list of indices
%%                into e of the associated edges (Nx3)
%%         et:    for each edge, a pair of indices of the triangles
%%                that are adjacent to this edge. One of the
%%                indices will be 0 for a boundary edge.
%%
%% Example: 
%%    fd=inline('sqrt(sum(p.^2,2))-1','p');
%%    [p,t]=distmesh2d(fd,@huniform,0.2,[-1,-1;1,1],[]);
%%    [e,te,et]=edgelist(p,t);
%%
%% DAM 11/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [e, te, et] = edgelist( p, t )

  N = size( t, 1 );
  te = zeros( N, 3 );

  % Long edge list contains for each edge
  % [ p0 p1 f e] where p0 < p1 are the indices of the endpoints, f is the
  % index in t that this edge came from, and e is an index labeling which
  % of the three edges this is.  For example, if [ 3 9 7 ] is the 43rd 
  % entry in t, we will generate three entries in the long edge list:
  %
  % [ 3 9 43 1 ]
  % [ 7 9 43 2 ]
  % [ 3 7 43 3 ]
  long_e = zeros( N*3, 4);

  % fill in long edge list
  for i=1:N 
    for j=1:3 
      M = t(i,j);
      m = t(i,mod( j, 3)+1);
      if( M < m )
        tmp = M;
        M = m;
        m = tmp;
      end
      eindex = 3*(i-1)+j; 
      long_e(eindex, 1) = m;
      long_e(eindex, 2) = M;
      long_e(eindex, 3) = i;
      long_e(eindex, 4) = j;
    end
  end

  % Sort the long edge list.
  long_e = sortrows( long_e, [1 2] );
  % Consolidate duplicate entries in the long edge list.  Before we begin
  % there will be two copies of each interior edge.
  k = 1;
  te( long_e(1,3), long_e(1,4) ) = 1;
  long_e( 1, 4 ) = 0;
  for i=2:3*N
      if( (long_e(i,1) == long_e(i-1,1)) & (long_e(i,2)==long_e(i-1,2))  )
        te( long_e(i,3), long_e(i,4) ) = k;
		long_e( k, 4 ) = long_e(i, 3 );
      else
        k=k+1;
        te( long_e(i,3), long_e(i,4) ) = k;
		long_e(k,:) = long_e(i,:);
		long_e(k,4) = 0;
      end
	  i = i + 1;
  end
  e = long_e(1:k,1:2);
  et = long_e(1:k,3:4);
  