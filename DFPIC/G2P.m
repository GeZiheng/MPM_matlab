function [vp,div_vp,Sp] = G2P(xi,vi,xp,kernel,transfer,index)
% transfer from grid to particles
[~,nx,ny] = size(xi);
[~,nPts] = size(xp);
vp = zeros(2,nPts);
div_vp = zeros(1,nPts);
xmin = xi(1,1,1);
ymin = xi(2,1,1);
h = xi(1,2,1)-xi(1,1,1);
nBss = length(index);
Sp = zeros(nBss,nPts);
for p = 1:nPts
    [bn_x,wx,dwx] = Kernel((xp(1,p)-xmin)/h,kernel);
    [bn_y,wy,dwy] = Kernel((xp(2,p)-ymin)/h,kernel);
    for i = 1:length(wx)
        for j = 1:length(wy)
            if bn_x+i>0 && bn_x+i<=nx && bn_y+j>0 && bn_y+j<=ny
                wip = wx(i)*wy(j);
                dwip = [dwx(i)*wy(j);wx(i)*dwy(j)]/h;
                vp(:,p) = vp(:,p) + wip * vi(:,bn_x+i,bn_y+j);
                % update matrix Sp
                Sp(:,p) = EvolveS(Sp(:,p),wip,vi(:,bn_x+i,bn_y+j),xi(:,bn_x+i,bn_y+j),xi(:,bn_x+2,bn_y+2),xp(:,p),h,transfer,index);
                % disp(Sp(:,p));
                % divergence of velocity
                div_vp(p) = div_vp(p) + vi(:,bn_x+i,bn_y+j)'*dwip;
            end
        end
    end
end