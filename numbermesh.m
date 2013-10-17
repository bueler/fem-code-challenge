function numbermesh(p,t,h0,unknown);
%NUMBERMESH  Adds numbering to nodes, triangles and unknown nodes in a mesh
%produced by distmesh2d and poissonv2 (or poissonDN).
%Example:
% >> fd=inline('sqrt(sum(p.^2,2))-1','p');  h0=0.3;
% >> [p,t]=distmesh2d(fd,@huniform,h0,[-1,-1;1,1],[]);
% >> numbermesh(p,t,h0)
%Example (continues above) including numbering of unknowns:
% >> f=inline('4','p');  [uh,in]=poissonv2(f,fd,h0,p,t);
% >> numbermesh(p,t,h0,in)
%   See also: POISSONV2, POISSONDN, DISTMESH2D, TRIMESH.
%ELB 11/13/04

clf, trimesh(t,p(:,1),p(:,2),zeros(size(p,1),1))
view(2), axis equal, axis off

for j=1:size(t,1)   % number triangles
    h=text(sum(p(t(j,:),1))/3,sum(p(t(j,:),2))/3,num2str(j));
    set(h,'Color','green'), end
%number all nodes (black numerals)
for j=1:size(p,1), text(p(j,1)+.1*h0,p(j,2),num2str(j)), end

if nargin==4  % if numbering of unknowns is desired
    un = p(unknown>0,:);  kn = p(~(unknown>0),:);  N=sum(unknown>0);
    % markers: * for unknown nodes; o for known
    hold on, plot(un(:,1),un(:,2),'*',kn(:,1),kn(:,2),'o')
    %number unknown nodes in red
    for j=1:N
        h=text(un(j,1)-.25*h0,un(j,2),num2str(j));
        set(h,'Color','red','FontSize',14,'FontWeight','bold'), end
end
