%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Function: apost.m
%%
%% Purpose:  Computes an a-posteriori estimate for the Dirichlet
%%           for the Poisson equation.
%%
%% Input:  uh:     approximate solution (nx1)
%%         f:     right-hand side for -\Delta u = f
%%         p:     matrix of points (nx2)
%%         t:     matrix of indices into p defining a triangulation (Nx3)
%%         
%%
%% Output: err:   Nx1 vector of squared element residuals
%%         faceh: Nx1 vector of the approximate diameter of each element
%% 
%% You can obtain L^2 and H^1 error estimates (up to a constant
%% depending on the trianglutation) with the following computations:
%%
%% L2 error: sqrt( sum( err.*faceh.^4 ) )
%% H1 error: sqrt( sum( err.*faceh.^2 ) )
%%
%% Example: 
%%
%%  fd=inline('sqrt(sum(p.^2,2))-1','p');
%%  [p,t] = distmesh2d( fd, @huniform, 0.15, [-1,-1;1,1],[]);
%%  f = inline('ones(size(p,1),1)','p');
%%  uh=poissonv3( f, p, t )
%%  [err, faceh]=apost( uh, f, p, t );
%%  l2errorest = sqrt( sum( err.*faceh.^4 ) )
%%
%% DAM 12/20
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [err, faceh]=apost( uh, f, p, t )

[ e te et ] = edgelist( p, t);
ff = feval(f,p);

% compute integration on reference triangle
% I_{ab} = \int \phi_a \phi_b
I = (ones(3,3) + eye(3))/24;

% gradient vectors of basis functions on reference triangle
% X_{ib} = \partial_i \phi_b
X = [ [ -1 -1 ]' eye(2) ];

% Exterior normal vectors;
n= [ 0 1 -1  ; -1 1 0]';
% Vectors parallel to each side, used for computing boundary integrals.
V = [ 1 -1 0  ; 0 1 1]';

err = zeros(size(t,1),1);

% We have to wait until both halves of the contribution to the
% error from boundary integrals have been computed before we
% know the whole contribution from each edge.  So we store the
% intermediate computations in edge_err.
edge_err = zeros(size(e,1),1);

% We return the approximate diameter of each edge in faceh.  It
% is, in fact, the square root of the area, which is correct to
% a constant for a regular triangulation.
faceh = zeros(size(t,1),1);

for( j = 1:size(t,1) )
    T = t(j,:);
    J = [ p(T(2),:) - p(T(1),:) ; p(T(3),:)-p(T(1),:)]';
    g = J'*J;
    invg = inv(g);
    for(k=1:3)
        N = invg * n(k,:)' / sqrt(n(k,:)*invg*n(k,:)');
        edge_err(te(j,k)) = edge_err(te(j,k)) + uh(T(:))'*X'*N/2;
    end

    % Compute the contribution from the interior of the element.
    dV  = sqrt(abs(det(g)));
    faceh(j) = sqrt(dV);
    tf = [ ff(T(1)) ff(T(2)) ff(T(3)) ];
    err(j) = tf * I * tf' * dV;
end

edge_err(:,1) = edge_err.^2;
for j=1:size(e,1)
    % Only have contributions from interior edges.
    if( et(j,2)>0 )
        f1 = et(j,1); err(f1) = err(f1) + edge_err(j,1);
        f2 = et(j,2); err(f2) = err(f2) + edge_err(j,1);
    end
end
