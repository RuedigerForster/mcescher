function [CornerIndex]=CornerFinder(x,y)
[a,b]=size(x);
CornerIndex=zeros(1,b);
m=0;
for i=2:(b-1)
    if (y(1,i)-y(1,i-1))>[[y(1,i+1)-y(1,i-1)]/[x(1,i+1)-x(1,i-1)]]*(x(1,i)-x(1,i-1))        
     m=m+1;  
     CornerIndex(1,m)=i;    
    end
end
CornerIndex=CornerIndex(1,1:m);

    
    