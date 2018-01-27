function c = CalcFuncC(r,h)
% calculate function C(r) for cohesion
if r>h/2 && r<=h
    c = 32/(pi*h^9) * (h-r)^3 * r^3;
else if r>0 && r<=h/2
        c = 32/(pi*h^9) * ( 2*(h-r)^3 * r^3 - h^6/64 );
    else
        c = 0;
    end
end
