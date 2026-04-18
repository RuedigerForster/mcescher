function Y = SmoothEB(X,a)
% There are 20 lines core codes in ATEB algorithm.
% This innovative idea comes from Dave Have, US;
% Modified and Change to MATLAB language by Xinbo Liu in Central South University, China.
%
%Article
%Baseline correction of high resolution spectral profile data based on exponential smoothing
%Xinbo Liu; Zhimin Zhang; Yizeng Liang; Pedro F. M. Sousa; Yonghuan Yun; LingYu
%https://doi.org/10.1016/j.chemolab.2014.09.018
%https://www.sciencedirect.com/science/article/pii/S016974391400207X?via%3Dihub#bb0205
  
b = 1.0-a; sx= 1.0; sy = a;
   yi= 0; yi = sy*yi+sx*X(1,1);
   Y(1,1) =yi;
   n = length(X);
   for i = 2:1:n
     yi= a*yi+b*X(1,i);
     Y(1,i) = yi;	
   end;
   sx = sx/(1.0+a); 
   sy = sy/(1.0+a);
   yi =sy*yi+sx*X(1,n);
   for i = n-1:-1:1
     yi = a*yi+b*Y(1,i);
     Y(1,n) =yi;
     Y(1,i) = yi;
   end;
end
