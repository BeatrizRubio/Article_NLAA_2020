function w_n_i = func_pesos(a, b)

n = size(a,2);
w_n_i = zeros(n,n+1);
w_n_i(1,1) = a(1);
w_n_i(1,2) = b(1);
for n = 2:size(w_n_i,1)
  for i = 1:size(w_n_i,2)
      if i <= n
        w_n_i(n,i) =  w_n_i(n,i)+ (a(n)*(n-i+1)/(n)*w_n_i(n-1,i));
      end
      if i > 1
        w_n_i(n,i) =  w_n_i(n,i) + (b(n)*(i-1))/(n)*w_n_i(n-1,i-1);
      end
  end
end

