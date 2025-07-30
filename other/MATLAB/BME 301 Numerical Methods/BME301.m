n=10;
a=zeros(n+1,n+1);
for r = 0:n
    for c = 0:n
        if ((r+c)==0) 
            break; 
        end
        a(r+1,c+1) = 1/(r+c);
    end
    if (r==0) 
        continue; 
    end
    a(r+4,r+4)=10;
end