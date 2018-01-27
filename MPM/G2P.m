function [xp,vp,grad_vp,div_vp] = G2P(xi,vi,xp,vp,vin,flip,dt,mode,kernel)
% transfer from grid to particles (APIC method)
[~,nx,ny] = size(xi);
nPts = length(xp(1,:));
grad_vp = zeros(2,2,nPts);
div_vp = zeros(1,nPts);
xmin = xi(1,1,1);
ymin = xi(2,1,1);
h = xi(1,2,1)-xi(1,1,1);
for p = 1:nPts
    newx = [0;0];
    vpic = [0;0];
    vflip = vp(:,p);
    [bn_x,wx,dwx] = Kernel((xp(1,p)-xmin)/h,kernel);
    [bn_y,wy,dwy] = Kernel((xp(2,p)-ymin)/h,kernel);
    for i = 1:length(wx)
        for j = 1:length(wy)
            if bn_x+i>0 && bn_x+i<=nx && bn_y+j>0 && bn_y+j<=ny
                wip = wx(i)*wy(j);
                dwip = [dwx(i)*wy(j);wx(i)*dwy(j)]/h;
                % transfer velocity
                vpic = vpic + wip * vi(:,bn_x+i,bn_y+j);
                %vflip = vflip + wip * (vi(:,bn_x+i,bn_y+j)-vin(:,bn_x+i,bn_y+j));
                % update matrix Bp (for APIC)
                grad_vp(:,:,p) = grad_vp(:,:,p) + wip * vi(:,bn_x+i,bn_y+j) * (xi(:,bn_x+i,bn_y+j)-xp(:,p))';
                % particle advection
                % newx = newx + wip * (xi(:,bn_x+i,bn_y+j)+vi(:,bn_x+i,bn_y+j)*dt);      % FLIP advection
                div_vp(p) = div_vp(p) + dwip' * vi(:,bn_x+i, bn_y+j);
            end
        end
    end
    if strcmp(mode,'APIC')
        switch kernel
            case 1
                grad_vp(:,:,p) = grad_vp(p);
            case 2
                grad_vp(:,:,p) = 4/h^2 * grad_vp(:,:,p);
            case 3
                grad_vp(:,:,p) = 3/h^2 * grad_vp(:,:,p);
        end
        %Cp(:,:,p) = Cp(:,:,p) - 1/2 * trace(Cp(:,:,p)) * eye(2);
        xp(:,p) = xp(:,p) + dt * vpic;         % APIC advection
        vp(:,p) = vpic;
    end
    if strcmp(mode,'FLIP')
        xp(:,p) = newx;
        vp(:,p) = (1-flip) * vpic + flip * vflip;
    end
end