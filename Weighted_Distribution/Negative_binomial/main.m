clear all
%Experiments results presented in
%E. Mainar, J.M. Pe\~na, B. Rubio, 
%Accurate bidiagonal decomposition of collocation matrices of 
%weighted $\varphi$ - transformed systems (2020),
%Numerical Linear Algebra Appl. e2295. 

format long E

f = @(x) x;
g = @(x) 1-x;
n=10;
B=zeros(n);

t = zeros(1,n);
b = zeros(1,n);
for i=1:n
    t(i) = i/(n+1);
    b(i) = (-1)^i * randi(100);
end

%Collocation matrix of the negative binomial basis
for i=1:n
	for j=1:n
		B(i,j)= (1-t(i))*nchoosek(n-1,j-1)*  (f(t(i)))^(j-1) *(g(t(i)))^(n-j);
	end 
end

%Bidiagonal decomposition of the collocation matrix
BDB = negative_binomial(t,n,f,g); 

%Solve linear system Ax=b
SolBDB=transpose(TNSolve(BDB,b))
SolMB=B\transpose(b) 



%function TNSolve(B,b)
%Solves a TN linear system Ax=b, where B=BD(A). (see TNSolve of Plamen Koev https://math.mit.edu/~plamen/software/TNTool.html)

%Using this bidiagonal decomposition, we can also obtaine the inverse, eigenvalues and singular values using the
%functions presented in  https://math.mit.edu/~plamen/software/TNTool.html.


