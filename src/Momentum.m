function [P,L] = Momentum(mp,xp,vp,xc,Sp,h,i1,i2)
% compute angular momentum
P = [0;0];
L = 0;
nPts = length(mp);
b = h^2/4;
for p=1:nPts
    P = P + mp(p) * vp(:,p);
    L = L + mp(p) * ( (xp(1,p)-xc(1)) * vp(2,p) - (xp(2,p)-xc(2)) * vp(1,p) + b * (Sp(i2,p)-Sp(i1,p)) );
end