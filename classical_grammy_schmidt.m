function [ Q, R ] = classical_grammy_schmidt( A )
%CLASSICAL_GRAMMY_SCHMIDT The classical grammy schmidt factorization
%   A naive QR factorization using successive orthogonalization. The result
%   of this function is Q and R such that A == Q*R

[n, m] = size(A);
V = zeros(n);
Q = zeros(n);
R = zeros(n);

for j=1:n
   V(j,:) = A(j,:);
   for i=1:j-1
       R(i,j) = Q(i,:)*A(j,:)';
       V(j,:) = V(j,:) - R(i, j)*Q(i,:);
   end
   R(j,j) = -norm(V(j,:), 2);
   Q(j,:) = V(j,:)/R(j,j);
end

end

