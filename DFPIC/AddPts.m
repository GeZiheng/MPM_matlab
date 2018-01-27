function [xp,vp,mp,nPts] = AddPts(xc,r,vc,wc,rho,shape,h)
% initialize particles
nPts = 0;
xp = zeros(2,1);
vp = zeros(2,1);
sample = poissonDisc([512,512],h*60/r);      % sample particles
sample = sample * 2*r/512 - r;
for p = 1:length(sample(:,1));
    % check if the point is inside object
    x = sample(p,1);
    y = sample(p,2);
    if IsObject(x,y,r,shape);
        nPts = nPts + 1;
        xp(1,nPts) = xc(1) + x;
        xp(2,nPts) = xc(2) + y;
        vp(1,nPts) = vc(1) - wc * y;
        vp(2,nPts) = vc(2) + wc * x;
    end
end
mp = ones(1,nPts) * rho * (h^2/10);  % mass of particles