function fi = GridForces(nPts,xi,xp,mp,Fe,Fp,J,Vp,mu,lambda,ksi,gamma,kernel,texture,model,beta)
% compute grid forces
[~,nx,ny] = size(xi);
xmin = xi(1,1,1);
ymin = xi(2,1,1);
h = xi(1,2,1)-xi(1,1,1);
fi = zeros(2,nx,ny);
fp = zeros(2,sum(nPts));
sigma = zeros(2,2,sum(nPts));
tot_n = 0;
%% surface tension and adhesion
for n = 1:length(nPts)
    if strcmp(texture(n,:),'Water')
        [ni,np] = Normal(xi,xp,kernel);
        switch model
            case 'CSF'
                for k = 1:nx
                    for l = 1:ny
                        cur = 0;
                        [bn_x,wx,dwx] = Kernel((xi(1,k,l)-xmin)/h,kernel);
                        [bn_y,wy,dwy] = Kernel((xi(2,k,l)-ymin)/h,kernel);
                        for i = 1:length(wx)
                            for j = 1:length(wy)
                                if bn_x+i>0 && bn_x+i<=nx && bn_y+j>0 && bn_y+j<=ny
                                    dwip = [dwx(i)*wy(j);wx(i)*dwy(j)]/h;
                                    % compute curvature
                                    if norm(ni(:,bn_x+i,bn_y+j))>0
                                        cur = cur + dwip' * ni(:,bn_x+i,bn_y+j)/norm(ni(:,bn_x+i,bn_y+j));
                                        %cur = cur + dwip' * ni(:,bn_x+i,bn_y+j);
                                    end
                                end
                            end
                        end
                        fi(:,k,l) = fi(:,k,l) - gamma(n)* cur * ni(:,k,l) ;
                    end
                end
            case 'CSS'
                for p = tot_n+1:tot_n+nPts(n)
                    if norm(np(:,p)) > 0
                        sigma(:,:,p) = gamma(n) * ( norm(np(:,p))^2/2*eye(2) - np(:,p)*np(:,p)' );
                    end
                end
            case 'SPH'
                id2 = tot_n+1:tot_n+nPts(n);
                fp(:,id2) = Cohesion(xp(:,id2),mp(:,id2),np(:,id2),2*h,gamma(n));
                for m = 1:n-1
                    id1 = sum(nPts(1:m-1))+1:sum(nPts(1:m));
                    [fp1,fp2] = Adhesion(xp(:,id1),xp(:,id2),mp(id1),mp(id2),2*h,beta);
                    fp(:,id1) = fp(:,id1) + fp1;
                    fp(:,id2) = fp(:,id2) + fp2;
                end
        end
    end
end
%% elasticity
for n = 1:length(nPts)
    for p = tot_n+1:tot_n+nPts(n)
        % compute stress tensor
        switch texture(n,:)
            case 'Jello'
                FE = Fe(:,:,p);         % elastic deformation
                % polar decomposition
                [U,S,V] = svd(FE);
                R = U*V';
                JE = det(FE);
                P = 2*mu(n) * (FE-R) + lambda(n) * JE*(JE-1)*inv(FE');
                tmp = Vp(p) * P * FE';
            case 'Snow'
                FE = Fe(:,:,p);         % elastic deformation
                FP = Fp(:,:,p);         % plastic deformation
                % polar decomposition
                [U,S,V] = svd(FE);
                R = U*V';
                JE = det(FE);
                JP = det(FP);
                P = 2*mu(n)*exp(ksi(n)*(1-JP)) * (FE-R) + lambda(n)*exp(ksi(n)*(1-JP)) * JE*(JE-1)*inv(FE');
                tmp = Vp(p) * P * FE';
            case 'Water'
                sigma(:,:,p) = sigma(:,:,p) + mu(n) * (1-J(p)^(-lambda(n))) * eye(2);
                %D = (grad_vp(:,:,k)+grad_vp(:,:,k)')/2;
                %sigma = sigma + 2*ksi * D - 2/3*ksi * div_vp(k)*eye(2);
                tmp = Vp(p) * J(p) * sigma(:,:,p);
            case 'Sand'
                FE = Fe(:,:,p);         % elastic deformation
                FP = Fp(:,:,p);         % plastic deformation
                % polar decomposition
                [U,S,V] = svd(FE);
                eps = log(S);
                P = U * (2*mu(n)*inv(S)*eps + lambda(n)*trace(eps)*inv(S)) * V';
                tmp = Vp(p) * P * FE';
        end
        [bn_x,wx,dwx] = Kernel((xp(1,p)-xmin)/h,kernel);
        [bn_y,wy,dwy] = Kernel((xp(2,p)-ymin)/h,kernel);
        for i = 1:length(wx)
            for j = 1:length(wy)
                if bn_x+i>0 && bn_x+i<=nx && bn_y+j>0 && bn_y+j<=ny
                    wip = wx(i)*wy(j);
                    dwip = [dwx(i)*wy(j);wx(i)*dwy(j)]/h;
                    % compute force
                    fi(:,bn_x+i,bn_y+j) = fi(:,bn_x+i,bn_y+j) - tmp * dwip + wip * fp(:,p);
                end
            end
        end
    end
    tot_n = p;
end