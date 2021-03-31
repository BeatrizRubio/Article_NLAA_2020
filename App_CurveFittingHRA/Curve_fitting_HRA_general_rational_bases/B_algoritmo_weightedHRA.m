function [x,y] = B_algoritmo_weightedHRA(w,P,f,g,t)

% Algorithm of evaluation and subdivision presented in Section 4 of
 %E. Mainar, J.M. Pe\~na, B. Rubio, Evaluation and subdivision algorithms for general 
%classes of totally positive rational bases, Computer Aided Geometric Design, 
%Volume 81, 2020, 101900, ISNN 0167-8396.



m=length(P);
allt=t;

hold on;
set(gca, 'XColor', 'none', 'YColor', 'none')


Q=P;
d=w;
curva=zeros(size(allt, 2), 2);
s=1;

%B-algorithm curve
for t = allt    
    id = m+1;
    for i=1:m-1
        for j=1:m-i
            d(i+1,j)=g(t)*d(i,j)+f(t)*d(i,j+1);
            resp = id-(m+1-i);
            Q(id,1:2)=(g(t)*d(i,j)/d(i+1,j))*Q(resp,1:2)+(f(t)*d(i,j+1)/d(i+1,j))*Q(resp+1,1:2);
            id = id+1;
        end
    end
    curva(s,1) = Q(id-1, 1);
    curva(s,2) = Q(id-1, 2);
    x(s) = curva(s,1); y(s) = curva(s,2);
    s=s+1;
end
plot(x,y,'o')

end 

  