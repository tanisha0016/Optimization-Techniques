clc
clear all

% phase 1 (define parameters)
z=@(x1,x2) 3*x1-5*x2;
c1=@(x1,x2) x1+x2-6;
c2=@(x1,x2) x1+x2-4;
c=[3 -5]
A=[1 1;1 1];
b=[6;4];
n=size(A,1);
m=size(A,2);


% phase 2 plotting
x1 = 0:0.1:max(b ./ A(:,1));
figure;

for i = 1:n
    x22 = (b(i) - A(i,1) .* x1) ./ A(i,2);
    x22 = max(0, x22);
    plot(x1, x22);
    hold on;
end

%phase 3 intersection
pt=[]
A=[A;eye(2)];
b=[b;zeros(2,1)];
for i=1:size(A,1)
    for j=i+1:size(A,1)
        AA=[A(i,:);A(j,:)];
        bb=[b(i);b(j)];
        if(det(AA)~=0)
            X=inv(AA)*bb;
            if(X>=0)
                pt=[pt X];
            end
        end
    end
end
disp(pt)

%phase 4 all feasible points
FP=[]
for i=1:length(pt)
    pt1=pt(1,i);
    pt2=pt(2,i);
    if(c1(pt1,pt2)>=0 && c2(pt1,pt2)<=0)
        FP=[FP pt(:,i)]
        plot(pt1,pt2,'*r',"markerSize",10)     
    end
end
disp(FP)


%phase 5 optimal solution
if(isempty(FP))
            disp("no feasible solution")
else
    Z=[]
    for i=1:length(FP)
        cost=z(FP(1,i),FP(2,i))
        Z=[Z cost]
    end
    [optimal_value,index]=min(Z)
    optimal_sol=FP(:,index)
    hold off
end
