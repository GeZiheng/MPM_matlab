function [P,L] = Momentum(mp,xp,vp,Cp,Dp)
% compute angular momentum
P = [0;0];
L = 0;
nPts = length(mp);
for k=1:nPts
    P = P + mp(k) * vp(:,k);
    L = L + mp(k) * ( xp(1,k)*vp(2,k) - xp(2,k)*vp(1,k) + Dp*(Cp(2,1,k)-Cp(1,2,k)) );
end