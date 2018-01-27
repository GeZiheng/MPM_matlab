function [mi,vi] = P2G(xp,mp,vp,Sp,xi,kernel,transfer,index)
% transfer from particles to the grid
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
                % update grid momentum
                %mvi(:,bn_x+i,bn_y+j) = mvi(:,bn_x+i,bn_y+j) + wip * mp(p) * vp(:,p);    % PIC transfer
                tmp_v = GetVFromS(vp(:,p),Sp(:,p),xi(:,bn_x+i,bn_y+j),xi(:,bn_x+2,bn_y+2),xp(:,p),h,transfer,index);
                mvi(:,bn_x+i,bn_y+j) = mvi(:,bn_x+i,bn_y+j) + wip * mp(p) * tmp_v;    % PolyPIC transfer
            end
        end
    end
end
%% compute grid velocity
for i = 1:nx
    for j = 1:ny
        if mi(i,j)>0
            vi(:,i,j) = mvi(:,i,j) / mi(i,j);
        end
    end
end