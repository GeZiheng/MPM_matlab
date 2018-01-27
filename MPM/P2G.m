function [mi,vi] = P2G(xp,mp,vp,grad_vp,xi,mode,kernel)
% transfer from particles to the grid (APIC method)
[~,nx,ny] = size(xi);
nPts = length(xp(1,:));
xmin = xi(1,1,1);
ymin = xi(2,1,1);
h = xi(1,2,1)-xi(1,1,1);
mi = zeros(nx,ny);
mvi = zeros(2,nx,ny);
vi = zeros(2,nx,ny);
%% transfer mass and momentum
for p = 1:nPts
    [bn_x,wx,dwx] = Kernel((xp(1,p)-xmin)/h,kernel);
    [bn_y,wy,dwy] = Kernel((xp(2,p)-ymin)/h,kernel);
    for i = 1:length(wx)
        for j = 1:length(wy)
            if bn_x+i>0 && bn_x+i<=nx && bn_y+j>0 && bn_y+j<=ny
                wip = wx(i)*wy(j);
                % update grid mass
                mi(bn_x+i,bn_y+j) = mi(bn_x+i,bn_y+j) + wip * mp(p);
                % update grid momentum ( with matrix Bp, and Dp = (h^2/4) * I )
                if strcmp(mode,'APIC')
                    mvi(:,bn_x+i,bn_y+j) = mvi(:,bn_x+i,bn_y+j) + wip * mp(p) * ( vp(:,p) + grad_vp(:,:,p) * (xi(:,bn_x+i,bn_y+j)-xp(:,p)));    % APIC transfer
                end
                if strcmp(mode,'FLIP')
                    mvi(:,bn_x+i,bn_y+j) = mvi(:,bn_x+i,bn_y+j) + wip * mp(p) * vp(:,p);    % FLIP transfer
                end
            end
        end
    end
end
%sum((sum(mvi,3)),2)
%% compute grid velocity
for i = 1:nx
    for j = 1:ny
        if mi(i,j)>0
            vi(:,i,j) = mvi(:,i,j) / mi(i,j);
        end
    end
end