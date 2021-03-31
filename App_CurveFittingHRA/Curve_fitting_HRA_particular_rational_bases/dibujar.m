function errores = dibujar(a, b, f, g, l, t, bx, by)
n = size(a,2);
format long E

%Computation of BDA with HRA
BDA = QRfuncBDA(a,b,f,g,t,n,l);

[Q,R]=TNQR(BDA);
   
R1 = zeros(n+1,n+1);
for i=1:n+1
    for j=1:n+1
        R1(i,j) = R(i,j);
    end
end
  
%Control points
dx=transpose(Q)*transpose(bx);

for i=1:n+1
    d1x(i)=dx(i);
end

for i=n+2:l
     d2x(i)=dx(i);
end

dy=transpose(Q)*transpose(by);

for i=1:n+1
    d1y(i)=dy(i);
end

for i=n+2:l
    d2y(i)=dy(i);
end

Solbx=transpose(TNSolve(R1,d1x)); %Coordinate x of control point with HRA
Solby=transpose(TNSolve(R1,d1y)); %Coordinate y of control point with HRA

P = [Solbx Solby];

plot(bx,by,'+');

%Plot curve
[xb,yb] = B_algoritmo_weightedHRA(a,b,P,f,g,t);

%HRA errors
vax= zeros(1,l);
vay= zeros(1,l);

for i=1:l
    vax(i)=abs(bx(i)-xb(i));
    vay(i)=abs(by(i)-yb(i));
end 

VAMx=(sum(vax)/l);
VAMy=(sum(vay)/l);
errores =[VAMx VAMy];


%function [Q,C]=TNQR(B);
%
% Given B=BD(A), computes Q and C=BD(R), where A=QR is a QR decomp. of A
%
% Written September 29, 2004. Added the computation of Q, June 9, 2007
% 
% Copyright (c) 2004 Plamen Koev. See COPYRIGHT.TXT for more details.

% function x=TNSolve(B,b)
% Solves a TN linear system Ax=b, where B=BD(A).
%
% This yields Bjorck-Pereyra type methods for every TN matrix A for which
% we have B=BD(A). A (and thus B) must be square
%
% Plamen Koev, Massachusetts Institute of Technology, December 23, 2005
%
% Copyright Plamen Koev, see COPYRIGHT.TXT for details.
