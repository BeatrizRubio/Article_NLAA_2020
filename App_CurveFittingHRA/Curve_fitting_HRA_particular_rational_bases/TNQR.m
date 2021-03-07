function [Q,B]=TNQR(B)

%function [Q,C]=TNQR(B);
%
% Given B=BD(A), computes Q and C=BD(R), where A=QR is a QR decomp. of A
%
% Written September 29, 2004. Added the computation of Q, June 9, 2007
% 
% Copyright (c) 2004 Plamen Koev. See COPYRIGHT.TXT for more details.

[m,n]=size(B);

Q=eye(m);

for i=1:n
   for j=m:-1:i+1
      x=B(j,i);
      B(j,i)=0;
      c=sqrt(1+x^2);
      B=TNAddToPrevious(B',x/c,c,j)';
      Q(:,[j-1,j])=Q(:,[j-1,j])*[1/c -x/c; x/c 1/c];
   end
end
