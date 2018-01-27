function fi = GridForces(xi,xp,J,Vp,mu,lambda,kernel)
% compute grid forces
[~,nx,ny] = size(xi);
[~,nPts] = size(xp);
xmin = xi(1,1,1);
ymin = xi(2,1,1);
h = xi(1,2,1)-xi(1,1,1);
fi = zeros(2,nx,ny);
for p = 1:nPts
    sigma = mu * (1-J(p)^(-lambda)) * eye(2);
    tmp = Vp(p) * J(p) * sigma;
    [bn_x,wx,dwx] = Kernel((xp(1,p)-xmin)/h,kernel);
    [bn_y,wy,dwy] = Kernel((xp(2,p)-ymin)/h,kernel);
    for i = 1:length(wx)
        for j = 1:length(wy)
            if bn_x+i>0 && bn_x+i<=nx && bn_y+j>0 && bn_y+j<=ny
                dwip = [dwx(i)*wy(j);wx(i)*dwy(j)]/h;
                % compute force
                fi(:,bn_x+i,bn_y+j) = fi(:,bn_x+i,bn_y+j) - tmp * dwip;
            end
        end
    end
end