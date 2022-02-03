function [x,C,P]=Gaussian_Elimination_4(A,b)
  
C=[A,b];  
[n]=size(C,1);
P=eye(n,n);

for j=1:(n-1)
    [p,q]=max(abs(C(j:n,j)));
    if q ~= 1 then
        C([j q+j-1],:)=C([q+j-1 j],:);
        P([j q+j-1],:)=P([q+j-1 j],:);
    end
    for i=(j+1):n
        C(i,j)=C(i,j)/C(j,j);
        C(i,j+1:n+1)=C(i,j+1:n+1)-C(i,j)*C(j,j+1:n+1);
    end
  end
  
x=zeros(n,1);

x(n)=C(n,n+1)/C(n,n);
for i=n-1:-1:1
    x(i)=(C(i,n+1)-C(i,i:n)*x(i:n))/C(i,i);
end

C=C(1:n,1:n);

endfunction
