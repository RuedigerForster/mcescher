LMV_RSA is a data analysis method used for chromatographic background drift correction

It can be performed in MATLAB environment

--------------------------------------------------------------------------------------------
Treat_GC.m and Treat_LCQTOF are two simple examples

--------------------------------------------------------------------------------------------

% function xbc = LMV_RSA(x,w)
%
% ------ Input Parameter ------
% x is the analyzed chromatogram
% w is the width of window, the default value is 30
%
% ------ Output Parameter ------
% xbc corrected chromatogram
