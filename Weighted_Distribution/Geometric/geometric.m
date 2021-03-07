function [BDA] = geometric(t,n)
% Bidiagonal decomposition of the collocation matrix of the Geometric basis
%E. Mainar, J.M. Pe\~na, B. Rubio, 
%Accurate bidiagonal decomposition of collocation matrices of 
%weighted $\varphi$ - transformed systems (2020),
%Numerical Linear Algebra Appl. e2295. 

format longE

BDA=zeros(n);

% Computation of the bidiagonal decomposition of the collocation matrix  

% Computation of the multipliers m_{i,j} of the collocation matrix of the
% Geometric basis
for i=2:n
    BDA(i,1)=t(i)/t(i-1);
end
 
for i=3:n
    for j=2:i-1
        aux=1;
        for k=1:j-1
            aux=(aux*(t(i-k)-t(i)))/(t(i-k-1)-t(i-1));
        end    
        BDA(i,j)=(t(i)/t(i-1))*aux;
    end
end


% Computation of the pivots p_{i,i} of the collocation matrix of the
% Geometric basis
BDA(1,1)=t(1);
for i=2:n
    aux=1;
    for k=1:i-1
        aux=aux*(t(k)-t(i));
    end
    BDA(i,i)=t(i)*aux;
end 


% Computation of the multipliers tilde m_{i,i} of the collocation matrix of the
% Geometric basis
for j=2:n
    for i=1:j-1   
        BDA(i,j)=(1-t(i));
    end
end
