%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%    File: findinterior.m
%%
%% Purpose: Returns a list of which edges in a triangulation are
%%          edges.
%%
%% Input:  p:     matrix of points (nx2)
%%         e:     list of edges (mx2) giving indices into p
%%         te:    triangle to edge lookup table (Nx3)
%%
%% Output: int_point:   (nx1) logical array; entry is true if the
%%                      corresponding point is a point in the
%%                      interior of the triangluation.
%%         int_edge:    (mx1) logical array; entry is true if the 
%%                      corresponding edge is an interior edge (so
%%                      at least one endpoint is an interior point.
%%
%% Given a triangluation specified by points p and triangles t, the
%% input arguments e and et can be obtained using 'edgelist.m'.
%
%% Example: 
%%
%%  fd=inline('sqrt(sum(p.^2,2))-1','p');
%%  [p,t] = distmesh2d( fd, @huniform, 0.15, [-1,-1;1,1],[]);
%%  [ e, te, et] = edgelist( p, t );
%%  [ ip, ie] = findinterior( p, e, te );
%%
%% DAM 11/24
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ int_point int_edge ] = findinterior( p, e, te )

% Find boundary edges by traversing edges and counting
bdy_count = zeros( size(e,1), 1 );
for( i=1:size(te,1) )
    for( j=1:3 )
        bdy_count( te(i,j) ) = bdy_count( te(i,j) ) + 1;
    end
end
% We traverse interior edges twice and boundary edges once.  All other
% edges are spurious.
bdy_edge = bdy_count == 1;
int_edge = bdy_count > 1;

% Find interior points as endpoints of interior edges that are not
% also endpoints of boundary edges.  We can't define the interior
% points as the complement of the boundary points since there might
% be extra points in the point list.
int_point=logical(zeros( size( p, 1 ), 1 ));
int_point( e(int_edge,1) ) = true;
int_point( e(int_edge,2) ) = true;
int_point( e(bdy_edge,1) ) = false;
int_point( e(bdy_edge,2) ) = false;
  