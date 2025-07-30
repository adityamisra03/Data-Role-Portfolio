function X=backsub(A,b)
%A is an n x n upper-triangular nonsingular matrix
%b is an n x 1 matrix
%X is the solution of Ax=b
n=max(size(b));
X=zeros(n,1);
X(n)=b(n)/A(n,n);
for k=n-1:-1:1
    X(k)=(b(k)-A(k,k+1:n)*X(k+1:n))/A(k,k);
end
