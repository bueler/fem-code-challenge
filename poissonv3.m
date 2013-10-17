%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Function: poissonv3.m
%%
%% Purpose:  Computes the an approximate solution to -\Delta u = f
%%           with u = 0 on the boundary of the domain.  Uses
%%           piecewise linear elements.  No display.
%%
%% Input:  f:     right-hand side for -\Delta u = f
%%         p:     matrix of points (nx2) that are vertices of a triangluation
%%         t:     matrix of indices into p defining a triangulation (Nx3)
%%         
%%
%% Output: uh:   The values of the approximate solution at each
%%                point p.
%%
%% Example: 
%%
%%  fd=inline('sqrt(sum(p.^2,2))-1','p');
%%  [p,t] = distmesh2d( fd, @huniform, 0.15, [-1,-1;1,1],[]);
%%  f = inline('ones(size(p,1),1)','p');
%%  uh=poissonv3( f, p, t );
%%
%% DAM 12/20
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uh=poissonv3( f, p, t )

[ e, te ] = edgelist( p, t);
ind = findinterior( p, e, te );
Np = size( p, 1); N = sum ( ind );
in = zeros( Np, 1 ); in(ind) = (1:N)';

A = sparse( N, N ); b = zeros( N, 1 );
ff = feval(f,p);

% compute integration on reference triangle
I = (ones(3,3) + eye(3))/12;

% gradient vectors of basis functions on reference triangle
X = [ [ -1 -1 ]' eye(2) ];

for( j = 1:size(t,1) )
    T = t(j,:);

    % Construct scaled metric on reference triangle.
    q(1,:) = p(T(2),:) - p(T(1),:);
    q(2,:) = p(T(3),:) - p(T(1),:);
    adT = [ q(2,2), -q(2,1) ; -q(1,2), q(1,1) ];
    d = abs(q(2,2)*q(1,1)-q(2,1)*q(1,2));
    g = adT * (adT') / d;

    i = in( T(1) ) ; j = in( T(2) ) ; k = in( T(3) );
    tf = [ ff(T(1)) ff(T(2)) ff(T(3)) ];
    C = g * X ;

    if i > 0
        A(i,i) =  A(i,i) + X(1,1)*C(1,1)+X(2,1)*C(2,1);
        b(i) = b(i) + tf * I(:,1) * d;
    end
    if j > 0
        A(j,j) =  A(j,j) + X(1,2)*C(1,2)+X(2,2)*C(2,2);
        b(j) = b(j) + tf * I(:,2) * d;
    end
    if  k > 0
        A(k,k) =  A(k,k) + X(1,3)*C(1,3)+X(2,3)*C(2,3);
        b(k) = b(k) + tf * I(:,3) * d;
    end

    if i*j > 0
        A(i,j) = A(i,j) + X(1,1)*C(1,2)+X(2,1)*C(2,2);
        A(j,i) = A(i,j);
    end
    if i*k > 0
        A(i,k) = A(i,k) + X(1,1)*C(1,3)+X(2,1)*C(2,3);
        A(k,i) = A(i,k);
    end
    if j*k > 0
        A(j,k) = A(j,k) + X(1,2)*C(1,3)+X(2,2)*C(2,3);
        A(k,j) = A(j,k);
    end
end

uh = zeros( Np , 1);
uh(ind) = A \ b;
