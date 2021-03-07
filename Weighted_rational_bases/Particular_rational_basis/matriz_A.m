function M =matriz_A(a,b,f,g,n,t)

omega=func_omega(a,b,f,g,n);
w=func_pesos(a,b);

M=zeros(n+1);

%Calculo matriz de colocación

for i=1:n+1
	for j=1:n+1
		M(i,j)= (w(n,j)*nchoosek(n,j-1)*(f(t(i)))^(j-1) *(g(t(i)))^(n+1-j))/omega(t(i));
     end   
end