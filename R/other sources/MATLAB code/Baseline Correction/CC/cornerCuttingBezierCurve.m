function [bezierCurve,FinalPoints,NewPoints]=cornerCuttingBezierCurve(x,y,steps)

[x2,y2]=CornerCutting(x,y,steps);

[FinalPoints,NewPoints]=AddMiddlePoints(x2,y2);

bezierCurve = bernsteinMatrix(size(FinalPoints,2)-1, linspace(0,1,5250))*FinalPoints';

% figure
% 
% plot(x,y,'b-')
% 
% hold on 
% 
% plot(bezierCurve(:,1),bezierCurve(:,2),'r-')