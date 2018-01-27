function cohesionF = Cohesion(xp,mp,np,d,gamma)
% calculate cohesion force
nPts = size(xp,2);
cohesionF = zeros(2,nPts);
for p = 1:nPts
    index = rangesearch(xp',xp(:,p)',d);
    index = cell2mat(index);
    neigh_num = size(index,2);
    if neigh_num == 1
        continue;
    end
    for k = 2:neigh_num
        q = index(1,k);
        dx = xp(:,p) - xp(:,q);
        cohesionF(:,p) = cohesionF(:,p) - gamma*mp(p)*mp(q) * CalcFuncC(norm(dx),d) * dx/norm(dx) + gamma*mp(p)*mp(q) * (np(:,p)-np(:,q));
    end
end