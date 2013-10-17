%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Function: adapt_poisson.m
%%
%% Purpose:  Computes the an approximate solution to -\Delta u = f
%%           with u = 0 on the boundary of the domain.  Uses
%%           piecewise linear elements and adaptive mesh refinement
%%           based on L^2 errors.
%%
%% Input:  f:     right-hand side for -\Delta u = f
%%         p:     matrix of points (nx2) that are vertices of a triangluation
%%         t:     matrix of indices into p defining a triangulation (Nx3)
%%         tol:   an L^2 error tolerance.  If 'u' is the true
%%                solution and uh is the approximate solution, then
%%                ||u-uh||_{L^2} < c tol where c is a constant
%%                of tol and f.  This last bit is lame; it would be
%%                ideal to figure out exactly what c is.              
%%
%% Output: uh:   The values of the approximate solution at each
%%                point p.
%%          p:   The adaptively refined mesh vertices 
%%          t:   The adaptively refined mesh elements
%%
%% Example: 
%%    p = [ [ 0 0 ]; [ 1 0 ]; [0 1] ; [1 1 ] ];
%%    t = [ [ 1 2 3 ]; [ 2 4 3 ] ];
%%    [ p t ] =  refine( p, t );
%%    f=inline('(1-sum(p.^2,2)).^5','p');
%%    [u,p,t]=adapt_poisson( f,p,t,0.0002);
%%
%% DAM 12/20
%% See also: POISSONV3, APOST, SREFINE, REFINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u,p,t]=adapt_poisson( f, p, t, tol )

tol2 = tol*tol;
while true
    u=poissonv3(f,p,t);
    [ err, faceh ]=apost( u, f, p, t );
    % Compute L2 error.
    err = err .* faceh.^4;
    vol=faceh.^2;
    totalVol=sum(vol,1);
    ploterr = err./vol*totalVol/tol2;

    % Trisurf is broken with regards to specifiying true color
    % patch colors.  So we inline some code here stolen
    % from trisurf that fixes the broken part.
    colors = [ ones(size(err)), 1-ploterr, 1-ploterr ];
    ax = newplot;
    h = patch('faces',t,'vertices',[p(:,1) p(:,2) zeros(size(p,1),1)],...
        'facevertexcdata', colors,...
        'facecolor',get(ax,'defaultsurfacefacecolor'), ...
        'edgecolor',get(ax,'defaultsurfaceedgecolor') );
    view(2),axis equal,axis off,drawnow

    % Change this next line to
    % if max(err./vol*totalVol) <= tol
    % to change the stopping condition to be uniform error density.
    if sum(err) <= tol2
        break;
    end
    [p,t]=srefine(p,t,err>tol2*vol/totalVol );
end
