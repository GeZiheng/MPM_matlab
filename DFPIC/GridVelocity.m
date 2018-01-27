function vi = GridVelocity(mi,vi,fi,di,g,dx,dt)
% update grid velocity
[nx,ny] = size(mi);
for i=1:nx
    for j=1:ny
        if mi(i,j) > 0
            % update grid velocity according to forces
            vi(:,i,j) = vi(:,i,j) + (fi(:,i,j)/mi(i,j)-[0;g]) * dt;
            % grid collision
            if di(i,j) < 2.2*dx
                n = -[di(i+(i<nx),j)-di(i+(i<nx)-1,j);di(i,j+(j<ny))-di(i,j+(j<ny)-1)];
                n = n / norm(n);
                vn = vi(:,i,j)' * n;
                if vn>0
                    vi(:,i,j) = vi(:,i,j) - vn * n;
                end
            end
        end
    end
end