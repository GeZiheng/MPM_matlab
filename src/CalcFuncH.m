function [l,dl] = CalcFuncH(w,p,h)
% calculate function H(r) for dfpic
l = w^2 + (2*p^3/h^2 - p/2)*w - h^2/4;
dl = 2*w + 2*p^3/h^2 - p/2;