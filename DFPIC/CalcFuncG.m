function [g,dg] = CalcFuncG(w,p,h)
% calculate function G(r) for polypic
g = w^2 + (4*p^3/h^2 - p)*w - h^2/4;
dg = 2*w + 4*p^3/h^2 - p;