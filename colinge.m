function [uh, p, t] = colinge;

% COLINGE FE implementation of Colinge and Rappaz (1999)
%   [uh, p, t] = colinge computes a finite element solution for a first order
%   approximation of the flow of a non-linear viscous fluid down an 
%   inclined slab. 
%
%   M. Truffer, Oct. 2004
%
%   see also: DISTMESH2D, FIXMESH

L=10;
H=2;                            % dimensions of the model
h0=.25;                         % average mesh size
T02=1e-10;                      % finite viscosity parameter
n=3;                            % flow law exponent
tol=1e-3;                       % convergence tolerance

% basal boundary condition

u_basal=inline(['0.5*exp(-(8*(x-' num2str(L) '/2)/' num2str(L) ').^2)']);

% mesh a rectangular domain using distmesh2d

figure(1)
fd=inline(['drectangle(p,0,' num2str(L), ',0,' num2str(H) ')']);
[p, t]= distmesh2d(fd,@huniform,h0,[0,0;L,H],[0,0;0,H;L,0;L,H]);

% distmesh2d can produce some double points, delete those

[p,t]=fixmesh(p,t);
Np=size(p,1); 
disp(['Mesh: ' num2str(Np) ' nodes; ' num2str(size(t,1)) ' elements'])

% find bdy vertices
geps=.001*h0;  
ind = (feval(fd,p) > -geps);  

% initial condition

u0=2/(n+1)*(1-(1-p(:,2)/2).^(n+1));

% Dirichlet bdy conditions

y=0:h0/100:H;                           % upstream bdy
Gam1(:,2)=y';
Gam1(:,1)=zeros(size(Gam1,1),1);
u1=2/(n+1)*(1-(1-y/2).^(n+1));

y=0:h0/100:H;                           % downstream bdy
Gam2(:,2)=y';
Gam2(:,1)=L*ones(size(Gam2,1),1);
u2=2/(n+1)*(1-(1-y/2).^(n+1));

x=0:h0/100:L;                           % basal bdy
Gam3(:,1)=x';
Gam3(:,2)=zeros(size(Gam3,1),1);
u3=u_basal(x);

% Neumann bdy cond.

Gam4(:,1)=(h0:h0/100:(L-h0))';          % top surface
Gam4(:,2)=ones(size(Gam4,1),1);

% go through all bdy points, assign the proper values

u_bound=zeros(1,Np);
for i=Np:-1:1
    if ind(i)
        for ii=1:4
            eval(['Gam=Gam' num2str(ii) ';']);
            [m(ii) mind(ii)]=min((Gam(:,1)-p(i,1)).^2+(Gam(:,2)-p(i,2)).^2);
        end
        [dummy, ii]=min(m);
        
        % points along Gam4 should be solved for (Neumann bdy cond)
        
        if ii~=4
            eval(['ub=u' num2str(ii) ';']);
            u_bound(i)=ub(mind(ii));
        else
            ind(i)=0;
        end
    end
end

go=1;

% Loop through all triangles

while go

    K=sparse(Np,Np);  F=zeros(Np,1);

    for q=1:size(t,1)
        % identify triangle vertices
        j=t(q,1); k=t(q,2); l=t(q,3);
    
        % following Reddy (1993, p. 434/435)
        % Each triangle is transformed to a "master element with corners
        % at (0,0) (1,0) and (0,1)
    
        b(1)=p(k,2)-p(l,2); b(2)=p(l,2)-p(j,2); b(3)=p(j,2)-p(k,2);
        c(1)=-p(k,1)+p(l,1); c(2)=-p(l,1)+p(j,1); c(3)=-p(j,1)+p(k,1);
        twoA=b(2)*c(3)-c(2)*b(3);           % Jacobian
        psix=b/twoA; psiy=c/twoA;
        ux=[u0(j) u0(k) u0(l)]*psix';       % gradient of u
        uy=[u0(j) u0(k) u0(l)]*psiy';
    
        % find the viscosity function by finding the zero of Eqn (16)
        % in Colinge and Rappaz
    
        s=ux^2+uy^2;
        v_guess=s^((n-1)/(2*n));
        visc=fzero(inline(['x^' num2str(2*n/(n-1)) '-' num2str(T02) '*x^2-' ...
                num2str(s)]),v_guess);
     
        % find the local stiffness matrix
    
        Kloc=(b'*b+c'*c)/(2*twoA*visc);
        floc=twoA/12*ones(3,1);         % load vector for f=1/2
    
        % assemble the global stiffness matrix
    
        K(t(q,:),t(q,:))=K(t(q,:),t(q,:))+Kloc;
        F(t(q,:))=F(t(q,:))+floc;
    end
    
    % remove Dirichlet bdy cond points: 
    % because Dirichlet bdy values are already known, the corresponding
    % rows and columns are deleted from the stiffness matrix. However, the
    % load vector needs to be adjusted accordingly.
    
    for i=Np:-1:1
        if ind(i)
            K(i,:)=[];
            F(i)=[];
            F=F-u_bound(i)*K(:,i);
            K(:,i)=[];            
        end
    end
 
    % solve for unknown u
    
    uh=K\F;
    
    % put the Dirichlet bdy pts back in
    
    for i=1:Np
        if ind(i)
            uh=[uh(1:i-1);u_bound(i);uh(i:size(uh,1))];
        end
    end
    
    % check tolerance
    
    if (norm(uh-u0))<tol
        go=0;
    else
        disp(['Convergence: ' num2str(norm(uh-u0))]);
        u0=uh;
    end
end

% plot results

x=linspace(0,L);
y=linspace(0,H);
[X Y]=meshgrid(x,y);
ZI = griddata(p(:,1),p(:,2),uh,X,Y);
surf(X,Y,ZI)
view(0,90)
colorbar('vert')
ch=get(gca,'Children');
set(ch,'LineStyle','none');

[err errind]=max(abs(uh-2/(n+1)*(1-(1-p(:,2)/2).^(n+1))));
disp(['Max error ' num2str(err) ' at ' num2str(p(errind,:))]);