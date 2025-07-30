function dy= FinalExam2a(x,y)

dy=zeros(2,1); %dy is a column vector!
dy(1) = y(2);
dy(2) = (R/L)*y(2)+(I/(L*C))*y(1)-(F0/L)*heaviside(t-1);
%end of function 
end