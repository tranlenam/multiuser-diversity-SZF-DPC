function [R,Q]=grams(A)
[m, n] = size(A);
R(1,1)=norm(A(1,:));
Q(1,:)=A(1,:)/R(1,1);
for j=2:m
    v=A(j,:);
    for k=1:j-1
        R(j,k)=A(j,:)*Q(k,:)';
        v=v-R(j,k)*Q(k,:);
    end
    if(norm(v)>(1e-10))
        R(j,j)=norm(v);
        Q(j,:)=v/norm(v);
%     else
%        break; 
    end
end 