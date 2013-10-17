%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%    File: srefine.m
%%
%% Purpose: Refines a subset of a  triangulation using edge bisection.
%%
%% Input:  p:     matrix of points (nx2)
%%         t:     matrix of indices into p defining a triangulation (Nx3)
%%         r:     vector of logicals indicating whcih triangles to
%%                refine (Nx1)
%%
%% Output: rp:    an augmented list of points in a refined triangluation.  
%%                The first N points will be the same as those in p.
%%         rt:    a list of triangles in the refined triangulation
%%
%% For example, see code in  apost. 
%%
%% DAM 11/24
%% See also: APOST, REFINE. 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rp, rt] = srefine( p, t, r )

	[e, te, et ] = edgelist( p, t );

	Np = size( p, 1);
	Nt = size( t, 1);
	Ne = size( e, 1);
	
	re = logical( zeros( Ne , 1 ) );

	% Determine what triangles to refine.  Triangles can have
    % either 1 or 3 edges bisected.  We require an interative loop
	% to make sure there are no triangles requiring two edges to be
	% bisected.
	while (true )
		ttype = zeros(Nt,1);
		re( te(r,1 ) ) = true;
		re( te(r,2 ) ) = true;
		re( te(r,3 ) ) = true;
		rei = find( re==true );
		
		for( j=1:size(rei))
			k=rei(j);
			ttype( et( k, 1) ) = ttype( et( k, 1) ) +1 ;
			if( et( k,2) ~=0 )
				ttype( et(k,2) ) = ttype( et(k,2) ) + 1; 
			end
		end

		% Make sure no 2's are in the ttype list
		if( sum( ttype==2, 1 ) == 0 )
			break
		end
		% If there are 2's, these triangles need to be completely refined.
		r = ttype > 1;
	end

	% Construct a list for looking up the refined edges.
	dp = sum( re,1 );
	rei     = zeros( Ne, 1 );
	rei(re) = 1:dp;

	% Construct the refined point list
	rp = zeros( size(p,1)+dp, 2 );
	rp( 1:size(p,1),:) = p;
	rp( Np+1:Np+dp,: ) = ( p( e(re,1),:) + p( e(re,2),:) ) /2;

	% The number of refined triangles can be computed from the ttype list
	Nrt = sum( (ttype + 1), 1 );
	rt = zeros( Nrt, 3 );
	
	% Copy the type 0 triangles over
	tloc = ttype==0;
	nt0 = sum(tloc,1);
	rt(1:nt0,:) = t(tloc,:);

	% Fully refine the type 3 triangles
	tloc = find( ttype == 3 );
	nt3 = size(tloc,1);
	
	for i=0:(nt3-1)
		k = tloc( i+1 );
		rpi = Np + rei( te( k, 1:3) );
		
		rt(nt0+4*i+1,:) = [ rpi(3), t(k,1), rpi(1) ];
		rt(nt0+4*i+2,:) = [ rpi(1), t(k,2), rpi(2) ];
		rt(nt0+4*i+3,:) = [ rpi(2), t(k,3), rpi(3) ];
		rt(nt0+4*i+4,:) = rpi;
	end
	
	% Partly refine the type 1 triangles
	tloc = find( ttype == 1 );
	nt1=size(tloc,1);
	
	for i=0:(nt1-1)
		k = tloc( i+1 );		
		refine_edge= find( re(te(k,:))==true );
		ind = mod( (-1:1)+refine_edge(1),3 ) + 1; 
		rpi = Np + rei( te( k, ind(1) ) );		
		rpp = t(k, ind);
		rt(nt0+4*nt3+2*i+1,:) = [ rpp(1), rpi, rpp(3) ];
		rt(nt0+4*nt3+2*i+2,:) = [ rpp(3), rpi, rpp(2) ];
	end
