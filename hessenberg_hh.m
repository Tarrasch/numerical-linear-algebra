function [ M ] = hessenberg_hh( A )
%HESSENBERG_HH Convert to hessenberg form using householder reflectors
%

[m, n] = size(A);
V = cell(1,m-2);

for k = 1:m-2
  x = A(k+1:m,k);
  v = [[sign(x(1)) * norm(x, 2)] zeros(1,m-k-1)]' + x;
  v = v/norm(v, 2);
  V{k} = v;
  A(k+1:m,k:m) = A(k+1:m,k:m) - 2*v*(v'*A(k+1:m,k:m));
  A(1:m,k+1:m) = A(1:m,k+1:m) - 2*(A(1:m,k+1:m)*v)*v';
end

M = A;
