clear all

%Experiments results presented in
%E. Mainar, J.M. Pe\~na, B. Rubio, 
%Accurate bidiagonal decomposition of collocation matrices of 
%weighted $\varphi$ - transformed systems (2020),
%Numerical Linear Algebra Appl. e2295. 

format long E


n=10;
A=zeros(n);

t = zeros(1,n);
b = zeros(1,n);

for i=1:n
    t(i) = i/(n+1);
    b(i) = (-1)^i * randi(100);
end

%Collocation matrix of the Poisson basis
for i=1:n
    for j=1:n
		A(i,j)=(t(i)^(j-1)/factorial(j-1))*exp(-t(i));
    end    
end

BDA = poisson(t,n); %Bidiagonal decomposition of the collocation matrix

%Solve linear system Ax=b

SolBDA=transpose(TNSolve(BDA,b))
SolMA=A\transpose(b)


%function TNSolve(B,b)
%Solves a TN linear system Ax=b, where B=BD(A). (see TNSolve of Plamen Koev https://math.mit.edu/~plamen/software/TNTool.html)

%Using this bidiagonal decomposition, we can also obtaine the inverse, eigenvalues and singular values using the
%functions presented in  https://math.mit.edu/~plamen/software/TNTool.html.

