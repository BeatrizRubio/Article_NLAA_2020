clear all
format longE


%Experiments results presented in
%E. Mainar, J.M. Pe\~na, B. Rubio, 
%Accurate bidiagonal decomposition of collocation matrices of 
%weighted $\varphi$ - transformed systems (2020),
%Numerical Linear Algebra Appl. e2295. 

%Particular rational basis. 

a=[1,1,1,1,1,1,1,1,1] %Change values of a. 
b=[1,3,9,2,2,2,1,1,1]

n = size(a,2);

t = zeros(1,n+1);
b = zeros(1,n+1);

%Choose functions f and g which satisfy the conditions for being Normalized B-basis. 
f = @(t) t;
g = @(t) 1 - t;

% f = @(t) t^2;
% g = @(t) 1 - t^2;

%Parameters of polynomial bases
for i=1:n+1
    t(i) = i/(n+2);
    b(i) = (-1)^i * randi(100);
end

% f = @(t) sin((pi/2-0.001+t)/2);
% g = @(t) sin((pi/2-0.001-t)/2);
% 
% f = @(t) sinh((pi/2-0.001+t)/2);
% g = @(t) sinh((pi/2-0.001-t)/2);

%Parameters of trigonometric and hyperbolic bases
% for i=1:n+1
%      t(i) = -pi/2 + 2*0.001 + i*pi/(n+1);
%      b(i) = (-1)^i * randi(100);
% end


%Bidiagonal decomposition of the collocation matrix. 
BDA=particular_rational(a,b,f,g,n,t) 
%Collocation matrix. 
M=matriz_A(a,b,f,g,n,t) 

%Solve linear system Ax=b
SolBD=transpose(TNSolve(BDA,b))
SolM = M\transpose(b)


%function TNSolve(B,b)
%Solves a TN linear system Ax=b, where B=BD(A). (see TNSolve of Plamen Koev https://math.mit.edu/~plamen/software/TNTool.html)

%Using this bidiagonal decomposition, we can also obtaine the inverse, eigenvalues and singular values using the
%functions presented in  https://math.mit.edu/~plamen/software/TNTool.html.


 
 



