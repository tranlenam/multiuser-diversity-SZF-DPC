function z=compvectors(x,y)
z=1;
n=max(length(x),length(y));
x=[x;zeros(n-length(x),1)];
y=[y;zeros(n-length(y),1)];
for k=1:n
    if(y(k)>x(k))
        z=0;
        break;
    end
end