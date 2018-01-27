function [xp,vp] = Advection(xp,vp,xi,dt)
% particle advection
xp = xp + vp * dt;
% particle collision
[~,nx,ny] = size(xi);
xmin = xi(1,1,1);
ymin = xi(2,1,1);
h = xi(1,2,1)-xi(1,1,1);
xmax = xmin + (nx-1) * h;
ymax = ymin + (ny-1) * h;
xp(1,:) = min(max(xp(1,:),xmin+1.8*h),xmax-1.8*h);
xp(2,:) = min(max(xp(2,:),ymin+1.8*h),ymax-1.8*h);
