function omega =func_omega(a,b,f,g,n)
w_n_i=func_pesos(a,b)
syms t;

for i = 1:n+1
    base(i) = w_n_i(n,i) * nchoosek(n,i-1) *f(t)^(i-1) * g(t)^(n-i+1);
end
omega=matlabFunction(sum(base));

