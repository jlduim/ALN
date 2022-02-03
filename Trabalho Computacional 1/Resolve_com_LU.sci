function [X]=Resolve_com_LU(C,B,P)
  
[n]=size(C,1);

X=[];
q=size(B,2);
for t=1:q
    h=(P')*B(:,t);
    x=zeros(n,1);
    y=zeros(n,1);
    y(1)=h(1);
    for i=2:n
        y(i)=(h(i)-C(i,1:i-1)*y(1:i-1));
    end
    x(n)=y(n)/C(n,n);
    for j=n-1:-1:1
        x(j)=(y(j)-C(j,j:n)*x(j:n))/C(j,j);
    end
    X=[X,x];
end

endfunction
