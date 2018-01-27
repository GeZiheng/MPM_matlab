function [ni,np] = Normal(xi,xp,kernel)
% compute normal vectors
[~,nx,ny] = size(xi);
xmin = xi(1,1,1);
ymin = xi(2,1,1);
h = xi(1,2,1)-xi(1,1,1);
ni = zeros(2,nx,ny);
nPts = length(xp(1,:));
np = zeros(2,nPts);
%% grid normal
for p = 1:nPts
    [bn_x,wx,dwx] = Kernel((xp(1,p)-xmin)/h,kernel);
    [bn_y,wy,dwy] = Kernel((xp(2,p)-ymin)/h,kernel);
    for i = 1:length(wx)
        for j = 1:length(wy)
            if bn_x+i>0 && bn_x+i<=nx && bn_y+j>0 && bn_y+j<=ny
                dwip = [dwx(i)*wy(j);wx(i)*dwy(j)]/h;
                ni(:,bn_x+i,bn_y+j) = ni(:,bn_x+i,bn_y+j) + dwip;
            end
        end
    end
end
%% particle normal
for p = 1:nPts
    [bn_x,wx,dwx] = Kernel((xp(1,p)-xmin)/h,kernel);
    [bn_y,wy,dwy] = Kernel((xp(2,p)-ymin)/h,kernel);
    for i = 1:length(wx)
        for j = 1:length(wy)
            if bn_x+i>0 && bn_x+i<=nx && bn_y+j>0 && bn_y+j<=ny
                wip = wx(i)*wy(j);
                np(:,p) = np(:,p) + wip * ni(:,bn_x+i,bn_y+j);
            end
        end
    end
end