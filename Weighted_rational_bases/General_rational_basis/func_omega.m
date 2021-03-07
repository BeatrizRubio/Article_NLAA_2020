function omega =func_omega(w,f,g,n)
syms t;

for i = 1:n+1
    base(i) = w(i) * nchoosek(n,i-1) *f(t)^(i-1) * g(t)^(n-i+1);
end
omega=matlabFunction(sum(base));

