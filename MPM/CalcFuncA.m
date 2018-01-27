function a = CalcFuncA(r,h)
% calculate function A(r) for adhesion
if r>h/2 && r<=h
    a = 0.7/h^3.25 * (-4*r^2/h+6*r-2*h)^(1/4);
else
    a = 0;
end
