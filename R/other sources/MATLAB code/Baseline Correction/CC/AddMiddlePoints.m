function [FinalPoints,NewPoints]=AddMiddlePoints(x2,y2)
[a,b]=size(x2);
NewPoints=zeros(2,b);
m=0;
k=1;
% plot(x2,y2,'b*-')
% 
% hold on 
for i=2:b-2  
    xnew=[x2(1,i)+x2(1,i+1)]/2;
    ynew=[y2(1,i)+y2(1,i+1)]/2;
    m=m+1;
    NewPoints(1,m)=xnew;
    NewPoints(2,m)=ynew;
    FinalPoints(1,k)=x2(1,i);
    FinalPoints(1,k+2)=x2(1,i+1);
    FinalPoints(2,k)=y2(1,i);
    FinalPoints(2,k+2)=y2(1,i+1);
        FinalPoints(1,k+1)=xnew;
    FinalPoints(2,k+1)=ynew;
    k=k+3;
%     plot(xnew,ynew,'r*')
    end


% plot(FinalPoints(1,:),FinalPoints(2,:),'b*-')
% 
% hold on 
% 
% plot(NewPoints(1,:),NewPoints(2,:),'r*')

%------------------------
% NewPoints=NewPoints(:,1:m);
% 
% FinalPoints=zeros(2,size(x2,2)+size(NewPoints,2));
% 
% FinalPoints(:,1:b)=[x2;y2];
% 
% FinalPoints(:,(b+1):end)=NewPoints;

% FinalPoints=FinalPoints';
% FinalPoints=sortrows(FinalPoints,1);
% 
% FinalPoints=FinalPoints';


    