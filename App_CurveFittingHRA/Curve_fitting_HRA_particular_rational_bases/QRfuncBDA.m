
function BDA = QRfuncBDA(a,b,f,g,t,n,l)

%In this algorithm we use the bidiagonal factorization for the squares collocation matrices
%of a particular class of rational function presented in  
%in Section 5 of 
%E. Mainar, J.M. Pe\~na, B. Rubio, 
%Accurate bidiagonal decomposition of collocation matrices of 
%weighted $\varphi$ - transformed systems (2020),
%Numerical Linear Algebra Appl. e2295. 

%We extend the mentioned bidigonal factorization to the case of rectangular 
%collocation matrices applying the least square method. 
    
omega=func_omega(a,b,f,g,n);
w=func_pesos(a,b);


BDA=zeros(l,n+1);


%Bidiagonal factorization
%Computation of the multipliers m_{i,j}


 for i=2:l
 	M=  g ( t( i ) )^(n)  / ( g( t(i-1) ) )^(n+1); 
	BDA( i,1)=M*g(t (i-1) ); 
    BDA( i,1)= BDA( i,1)*(omega(t(i-1))/omega(t(i)));  % cambio para racionales
	for j=2: i-1
      if (j<=n+1)
	     M=M* g( t ( i-1) )/  g( t(i) ) *  ( f(t(i) )* g(t( i- j+1) )  -  f( t( i-j+1 ) )*g( t( i) ) ) /  ( f( t( i-1) )*g( t(i- j) )-f( t( i-j) )*g( t(i-1) ) );   
	     BDA(i,j)= M * g( t (i-j) ); 
         BDA(i,j)= BDA(i,j)*(omega(t(i-1))/omega(t(i)));   % cambio para racionales
       end
	end 
 end

  %Computation of the pivots p_{i,i}
   BDA(1,1)=g(t(1))^(n);
   BDA(1,1)=BDA(1,1)*w(n,1)/(omega(t(1)));   
   q=1;
   for  i=2:n+1
  	q=q* (n-i+2)/(i-1)/g(t(i-1));
   	aux=1;
   	for k=1: i-1
   		aux=aux*   ( f( t(i) ) * g( t( k) ) - f( t(k) )*g( t(i) ) );
   	end 
        BDA(i,i)=q*aux*(g(t(i)))^(n-i+1);
        BDA(i,i)=BDA(i,i)*w(n,i)/(omega(t(i)));   
    end  
%Computation of the multipliers tilde m_{i,i}

  for j=1:n	
     coef= f(t(j))/g(t(j));
            for i= j+1:n+1 
    	BDA(j,i )=coef* ( n-i+2) / (i-1);
        BDA(j,i)=(w(n,i)/w(n,i-1))*BDA(j,i);     
    	end
  end

end
   
