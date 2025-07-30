function [Aug,X]=uptrbk(A,b)
%A is an N x N nonsingular matrix
%b is an N x 1 matrix
%X is an N x 1 matrix containing the solution of Ax=b
[~ , N]=size(A);
X=zeros(N,1);
C=zeros(1,N+1);
Aug=[A b];
for p=1:N-1
    [Y,j]=max(abs(Aug(p:N,p))); %notice Y is not used.  Why?
    C=Aug(p,:);
    Aug(p,:)=Aug(j+p-1,:); %swaps rows
    Aug(j+p-1,:)=C;
    if (Aug(p,p)==0) % cannot pivot
        disp('A singular, No solution');
    break;
    end
    for k=p+1:N
        m=Aug(k,p)/Aug(p,p);
        Aug(k,p:N+1)=Aug(k,p:N+1)-m*Aug(p,p:N+1);
    end
end
X=backsub(Aug(1:N,1:N),Aug(1:N,N+1));

end