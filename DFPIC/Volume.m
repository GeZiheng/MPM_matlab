function Vp = Volume(mi,xi,mp,xp,kernel)
% compute particle volumes
[d,nx,ny] = size(xi);
nPts = length(xp(1,:));
xmin = xi(1,1,1);
ymin = xi(2,1,1);
h = xi(1,2,1)-xi(1,1,1);
Vp = zeros(1,nPts);
for k = 1:nPts
    rho = 0;
    % compute density through interpolation
    [bn_x,wx,dwx] = Kernel((xp(1,k)-xmin)/h,kernel);
    [bn_y,wy,dwy] = Kernel((xp(2,k)-ymin)/h,kernel);
    for i = 1:length(wx)
        for j = 1:length(wy)
            if bn_x+i>0 && bn_x+i<=nx && bn_y+j>0 && bn_y+j<=ny
                wip = wx(i)*wy(j);
                rho = rho + mi(bn_x+i,bn_y+j) * wip / h^3;
            end
        end
    end
    Vp(k) = mp(k) / rho;
end   