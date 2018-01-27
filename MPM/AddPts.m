function [xp,vp,mp,nPts] = AddPts(c,r,v,w,rho,shape,h,mode)
% initialize particles
nPts = 0;
xp = zeros(2,1);
vp = zeros(2,1);
if strcmp(mode,'Poisson')
    sample = poissonDisc([512,512],h*60/r);      % sample particles
    sample = sample * 2*r/512 - r;
end
if strcmp(mode,'Uniform')
    [xx,yy] = meshgrid(-r:h/2.5:r,-r:h/2.5:r);
    sample = [reshape(xx,size(xx,1)*size(xx,2),1),reshape(yy,size(yy,1)*size(yy,2),1)];
end
for p = 1:length(sample(:,1));
    % check if the point is inside object
    x = sample(p,1);
    y = sample(p,2);
    if IsObject(x,y,r,shape);
        nPts = nPts + 1;
        xp(1,nPts) = c(1) + x;
        xp(2,nPts) = c(2) + y;
        vp(1,nPts) = v(1) + w * y;
        vp(2,nPts) = v(2) - w * x;
    end
end
mp = ones(1,nPts) * rho * (h^2/10);  % mass of particles