function [x2,y2]=CornerCutting(x,y,Steps)
x2=x;
y2=y;
for i=1:Steps    
    [CornerIndex]=CornerFinder(x2,y2);    
%     plot(x2,y2,'b*-')    
    x2(:,CornerIndex)=[];
    y2(:,CornerIndex)=[];  
end



