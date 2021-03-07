function [BDB] = negative_binomial(t,n,f,g)
% Bidiagonal decomposition of the collocation matrix of the negative
% binomial basis

%E. Mainar, J.M. Pe\~na, B. Rubio, 
%Accurate bidiagonal decomposition of collocation matrices of 
%weighted $\varphi$ - transformed systems (2020),
%Numerical Linear Algebra Appl. e2295. 
format longE

BDA=zeros(n);
BDB=zeros(n); 

% Computation of the bidiagonal decomposition of the collocation matrix  
% Firstly, computation of the bidiagonal decomposition of the collocation
% matrix of the Bernstein basis

% Computation of the multipliers m_{i,j} of the collocation matrix of the
% Berstein basis
for i=2:n 
	M= g ( t( i ) )^(n-1)  / ( g( t(i-1) ) )^(n); 
	BDA( i,1)=M*g(t (i-1) );
    for j=2:i-1
        M=M* g( t ( i-1) )/  g( t(i) ) *  ( f(t(i) )* g(t( i- j+1) )  -  f( t( i-j+1 ) )*g( t( i) ) ) /  ( f( t( i-1) )*g( t(i- j) )-f( t( i-j) )*g( t(i-1) ) );   
        BDA(i,j)= M * g( t (i-j) ) ;	   
    end
end

% Computation of the pivots p_{i,i} of the collocation matrix of the
% Berstein basis
BDA(1,1)=g(t(1))^(n-1);

q=1;
for i=2:n
    q=q*(n-i+1)/(i-1)/g(t(i-1));
    aux=1;
    for k=1: i-1
        aux=aux* ( f( t(i) ) * g( t( k) ) - f( t(k) )*g( t(i) ) );
    end 
    BDA(i,i)=q*aux*(g(t(i)))^(n-i);
end

% Computation of the multipliers tilde m_{i,i} of the collocation matrix of the
% Berstein basis
for j=1:n-1	
    coef= f(t(j))/g(t(j));
    for i= j+1:n 
        BDA(j,i)=coef*(n-i+1)/(i-1);
    end
end

% Finally, applying Theorem 4, computation of the bidiagonal decomposition
% of the collocation matrix of the negative binomial basis

% Computation of the pivots p_{i,i} of the collocation matrix of the
% negative binomial basis
for i=1:n
    BDB(i,i)=(1-t(i))*BDA(i,i);
end

% Computation of the multipliers tilde m_{i,i} of the collocation matrix of the
% negative binomial basis
for i=1:n
    for j=1:i-1
        BDB(i,j)=(1-t(i))/(1-t(i-1))*BDA(i,j);
    end
    for j=i+1:n
        BDB(i,j)=BDA(i,j);
    end
end