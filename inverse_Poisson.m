function [src,src_exact]=inverse_Poisson
%INVERSE_POISSON  Consider the following equation
%    -\Delta u = f, u\in[0,1]
%with boundary conditions u(0)=u(1)=0.  Suppose we know its solution u(x) 
%at the certain sets of points  {x_i}_{i=1}^n \in [0,1]  at which 
%u(x_i)=v_i.  We want to find the source function f(x), x \in [0,1].

%The idea is to minimize the functional |u(x_i)-v_i|^2, by A^*((Ag)(x_i)-v_i), 
%where A: (rhs of the above equation)->u(x_i), and * stands for the adjoint
%operator.

%space points
step=0.01;
x=0:step:1;
n=length(x);

%where we measure u(x)
v=[20 40 60 70 80];

k=n+2;

A=sparse(k,k);
F=zeros(k,1);
%the test source function
src_exact=x.*(1-x).*(0.5-x);
src=src_exact;

A(1,1)=-1/(x(2)-x(1));
A(1,2)= 1/(x(2)-x(1));
for i=2:n-1
    A(i,i-1)=1/(x(i)-x(i-1));
    A(i,i)=-1/(x(i)-x(i-1)) - 1/(x(i+1)-x(i));
    A(i,i+1)=1/(x(i+1)-x(i));
end
A(n,n)=-1/(x(n)-x(n-1));
A(n,n-1)=1/(x(n)-x(n-1));
A(k-1,n)=1;
A(n,k-1)=1;
A(k,1)=1;
A(1,k)=1;

F(1:n)=src*step;

F(k-1)=0;
F(k)=0;

X=A\F;

%got solutions, to determine src_exact
d=X(v);
src=zeros(n,1);

for iter=1:500
    %Direct problem
    F(1:n)=src*step;
    X=A\F;

    %Adjoint problem
    D=sparse(1,k);
    for vc=1:length(v)
        D(vc,v(vc))=1;
    end

    B=sparse(k,k);
    U=sparse(k,1);
    for m=1:k
        Em=zeros(k,1);
        Em(m)=1;
        B(m,:)=(A*Em)';
    
        U(m)=sum((D*Em).*(D*X-d));
    end

    Z=B\U;

    %Check the Lagrange identity
    lhs=sum(Z(1:n).*F(1:n));
    rhs=sum((D*X-d).*X(v));

    %Find the minimimum by Z_{n+1}=Z_n+\beta grad{Z_n}
    beta=sum((D*X-d).^2)/sum(Z(1:n).^2);
    src=src-beta*Z(1:n);
    
    %Compare exact and trial sources
    figure(1)
    plot(x,src,'-r'),  xlabel x,  ylabel('Source')
    hold on,  plot(x,src_exact,'-b'),  hold off
    title('red line is approximation, blue is exact')
    text(0.2,-0.03,num2str(norm([lhs rhs],'inf')))
    pause(0.1);
end

%figure(3),  plot(x,X(1:n)),  hold on,  plot(x,Z(1:n)),  hold off

