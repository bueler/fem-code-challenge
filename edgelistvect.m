%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%    File: edgelistvect.m
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
%%
%% Example: 
%%
%%  fd=inline('sqrt(sum(p.^2,2))-1','p');
%%  [p,t] = distmesh2d(fd, @huniform, 0.5, [-1,-1;1,1],[]);
%%  [e,te] = edgelistvect(p, t); 
%% 
%% DAM 10/26
%% Vectorized: Jed 2004-11-16
%% renamed: ELB 11/17/2004
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [e, te] = edgelistvect(p, t)

N = size(t,1);
te = zeros(N,3);

% Long edge list contains for each edge
% [ p0 p1 f e] where p0 < p1 are the indices of the endpoints, f is the
% index in t that this edge came from, and e is an index labeling which
% of the three edges this is.  For example, if [ 3 9 7 ] is the 43rd 
% entry in t, we will generate three entries in the long edge list:
%
% [ 3 9 43 1 ]
% [ 7 9 43 2 ]
% [ 3 7 43 3 ]
long_e = zeros(N*3,4);

% fill in long edge list using some index magic
i = reshape(repmat(1:N,[3 1]),[3*N 1]);
j = repmat((1:3)',[N 1]);
jm = repmat([2;3;1],[N 1]);
M = t(sub2ind(size(t),i,j));
m = t(sub2ind(size(t),i,jm));
long_e = [min(m,M) max(m,M) i j];

long_e = sortrows(long_e, [1 2]);                 % Sort the long edge list.
d = [1==1; sum(abs(diff(long_e(:,1:2))),2) ~= 0]; % Tag unique elements
e = long_e(d,1:2);                                % Extract unique elements
te(sub2ind(size(te), long_e(:,3), long_e(:,4))) = cumsum(d);
