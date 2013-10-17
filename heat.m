%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Function: heat.m
%%
%% Purpose:  Computes an approximate solution to the heat
%%           equation 
%%                        u_t-\Delta u = f
%%           with u=u_0 at time 0 and with u = 0 on the boundary of 
%%           the domain.  Uses piecewise linear elements and
%%           backwards Euler timestepping.
%%
%% Input: u0:     function giving the initial value of u
%%         f:     forcing function for the right hand side of the
%%                heat equation
%%         T:     last time solution should be computed up to
%%         n:     Number of time steps.
%%         p:     matrix of points (mx2) that are vertices of a triangluation
%%         t:     matrix of indices into p defining a triangulation(Nx3)
%%
%% Output: u:     An m x (n+1) matrix giving function values of the 
%%                approximate solution (m=#of points, n=#of timesteps)
%%
%% Example: 
%%
%% fd=inline('sqrt(sum(p.^2,2))-1','p');
%% [p,t] = distmesh2d( fd, @huniform, 0.1, [-1,-1;1,1],[]);
%% f=inline('8*p(:,1)*sin(t)+(1-sum(p.*p,2)).*p(:,1)*cos(t)','p','t');
%% u0=inline('zeros(size(p,1),1)','p');
%% uh = heat( u0, f, 2*pi, 100, p,t );
%%
%%
%% DAM 12/20
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = heat( u0, f, T, n, p, t )

	%% Compute the matrices A and B giving the inner producs in
    %% in the H^1_0 and L^2 norms.
    [ e, te ] = edgelist( p, t);
    ind = findinterior( p, e, te );
	Np = size( p, 1); N = sum ( ind );
	in = zeros( Np, 1 ); in(ind) = (1:N)';
	A = sparse( N, N ); B = sparse( N, N ); 

	% gradient vectors of basis functions on reference triangle
	X = [ [ -1 -1 ]' eye(2) ];
	% compute integration on reference triangle
	I = (ones(3,3) + eye(3))/24;

	for( j = 1:size(t,1) )
   		 Tr = t(j,:);
        % Construct scaled metric on reference triangle.
	    q(1,:) = p(Tr(2),:) - p(Tr(1),:); 
   	    q(2,:) = p(Tr(3),:) - p(Tr(1),:); 
		adT = [ q(2,2), -q(2,1) ; -q(1,2), q(1,1) ];
		d = abs(q(2,2)*q(1,1)-q(2,1)*q(1,2));
		g = adT * (adT') /2/d;

		i = in( Tr(1) ) ; j = in( Tr(2) ) ; k = in( Tr(3) );
		C = g * X ;

		if i > 0  
			A(i,i) =  A(i,i) + X(1,1)*C(1,1)+X(2,1)*C(2,1);
			B(i,i) =  B(i,i) + d/12;
    	end
		if j > 0  
           	A(j,j) =  A(j,j) + X(1,2)*C(1,2)+X(2,2)*C(2,2);
			B(j,j) = B(j,j) + d/12;
    	end
		if  k > 0 
        	A(k,k) =  A(k,k) + X(1,3)*C(1,3)+X(2,3)*C(2,3);
			B(k,k) = B(k,k) + d/12;
    	end

	    if i*j > 0 
    	    A(i,j) = A(i,j) + X(1,1)*C(1,2)+X(2,1)*C(2,2); A(j,i) = A(i,j);
			B(i,j) = B(i,j) + d/24; B(j,i) = B(i,j);
		end
    	if i*k > 0 
        	A(i,k) = A(i,k) + X(1,1)*C(1,3)+X(2,1)*C(2,3); A(k,i) = A(i,k); 
			B(i,k) = B(i,k) + d/24; B(k,i) = B(i,k);
	    end
    	if j*k > 0 
        	A(j,k) = A(j,k) + X(1,2)*C(1,3)+X(2,2)*C(2,3); A(k,j) = A(j,k);
			B(k,j) = B(k,j) + d/24; B(j,k) = B(k,j);
		end
	end		

	% Initialize the return value
	u=zeros( size(p,1),n+1);
	u(:,1) = feval(u0,p);

	% Construct a determinant vector
	q= p( t(:,2),: ) - p(t(:,1),:);
	r= (p( t(:,3),: ) - p(t(:,1),:)) * [ [ 0 -1 ] ; [ 1, 0 ] ];
	d= abs(sum( q.*r, 2));
	
	dt = T/n;	
	BA=B+dt*A;	

	% Time step.
	time=0;		
	for( l=1:n )
		F = zeros(sum(ind,1),1);			
		time = time + dt;
		ff = feval(f,p,time);
		for( m = 1:size(t,1) )
			i = in( t(m,1) ) ; j = in( t(m,2) ) ; k = in( t(m,3) );
			tf = [ ff(t(m,1)) ff(t(m,2)) ff(t(m,3)) ];
			if i > 0
				F(i) = F(i) + tf * I(:,1) * d(m);
			end
			if j > 0
				F(j) = F(j) + tf * I(:,2) * d(m);
			end
			if  k > 0
				F(k) = F(k) + tf * I(:,3) * d(m);
			end
		end
		u(ind,l+1)=BA \( B*u(ind,l)+dt*F );	
	end
