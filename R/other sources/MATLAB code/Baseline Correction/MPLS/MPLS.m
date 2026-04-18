function z= MPLS(x, lambda, windowWeight,order)
%  MPLS means Morphological Weighted Penalized Least Squares;		
%  Input 
%         x: vector with n variables
%         lambda: lambda is an adjustable parameter, it can be adjusted by user. 
%         The larger lambda is, the smoother z will be 
%         windowWeight:windowWeight is half the window of the structuring element 
%         order: an integer indicating the order of the difference of penalties
%         
%  Output
%         z: the fitted vector
%        
%  Reference:
%         (1) Eilers, P. H. C., A perfect smoother. Analytical Chemistry 75 (14), 3631 (2003).
%         (2) Eilers, P. H. C., Baseline Correction with Asymmetric Least
%         Squares Smoothing, http://www.science.uva.nl/~hboelens/publications/draftpub/Eilers_2005.pdf
%         (3) Zhimin. Zhang, Shan. Chen, and Yizeng. Liang, Baseline correction 
%             using adaptive iteratively reweighted penalized least squares. Analyst 135 (5), 1138-1146 (2010).stems 82 (1-2), 59 (2006).
%  Programmer: dejian zhan @ central south university on April 3,2013
if nargin ==3
    order=1;
end
if nargin ==2
    windowWeight=50;
end
Xbg=Morphology(x,windowWeight);
m = length(Xbg);
w = zeros(m, 1);
repeat = 1;
i=1;
Z=[];
xbg_reversel=Xbg(m:-1:1);
temp1=ones(m,1);
temp1(find([1;diff(Xbg, 1)]==0))=0;
temp2=ones(m,1);
temp2(find([1;diff(xbg_reversel, 1)]==0))=0;
temp2_reversel=temp2(m:-1:1);
temp=zeros(m,1);
temp3=temp1+temp2_reversel;
if temp3(2)==0
    temp3(1)=1;
end
temp(find(temp3==1))=1;
index=find(temp==1);
w = zeros(m, 1);
for i=1:floor(length(index)/2)-1
    [M,N]=min(Xbg(index(2*(i-1)+2):index(2*(i)+1)));
    w(index(2*(i-1)+2)+N-1)=1;
end
w(m)=5;
w(1)=5;
   z = WhittakerSmooth(Xbg, lambda, w, order);
   Z=[Z z];
end


function z = WhittakerSmooth(y, lambda, w, d)
% Weighted smoothing with a finite difference penalty
% y:      signal to be smoothed
% lambda: smoothing parameter 
% w:      weights (use one for missing values)
% d:      order of differences in penalty (generally 2)
% Paul Eilers, 2002
m = length(y);
W = spdiags(w, 0, m, m);
D = diff(speye(m), d);
C = chol(W+lambda * D' * D);
z = C \ (C' \ (w .* y));
end

function Xbg=Morphology(X,windowWeight)
if nargin ==1
    windowWeight=50;
end
Nx=length(X);
if windowWeight >=Nx/2-1
    windowWeight=floor(Nx/2-1);
end

Xerosion=zeros(1,Nx);
Xdilation=zeros(1,Nx);
i=1;

for j=1:windowWeight
    Xerosion(i,j)=min(X(i,1:j+windowWeight));
end
for j=windowWeight+1:(Nx-windowWeight)
    Xerosion(i,j)=min(X(i,j-windowWeight:j+windowWeight));
end
for j=(Nx-windowWeight+1):Nx
    Xerosion(i,j)=min(X(i,j-windowWeight:Nx));
end

for j=1:windowWeight
    Xdilation(i,j)=max(Xerosion(i,1:j+windowWeight));
end
for j=windowWeight+1:(Nx-windowWeight)
    Xdilation(i,j)=max(Xerosion(i,j-windowWeight:j+windowWeight));
end
for j=(Nx-windowWeight+1):Nx
    Xdilation(i,j)=max(Xerosion(i,j-windowWeight:Nx));
end
 Xbg=Xdilation';
 for j=1:windowWeight
    Xdilation(i,j)=max(X(i,1:j+windowWeight));
end
for j=windowWeight+1:(Nx-windowWeight)
    Xdilation(i,j)=max(X(i,j-windowWeight:j+windowWeight));
end
for j=(Nx-windowWeight+1):Nx
    Xdilation(i,j)=max(X(i,j-windowWeight:Nx));
end
Xbg(end-3:end)=X(end-3:end);
end