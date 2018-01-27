function [adhF1,adhF2] = Adhesion(xp1,xp2,mp1,mp2,d,beta)
% calculate adhesion force between xp1 and xp2
n1 = size(xp1,2);
n2 = size(xp2,2);
adhF1 = zeros(2,n1);
adhF2 = zeros(2,n2);
for p = 1:n1
    index = rangesearch(xp2',xp1(:,p)',d);
    index = cell2mat(index);
    neigh_num = size(index,2);
    if neigh_num == 0
        continue;
    end
    fprintf('bang!');
    for k = 1:neigh_num
        q = index(1,k);
        dx = xp1(:,p) - xp2(:,q);
        force = beta*mp1(p)*mp2(q) * calcFuncA(norm(dx),d) *dx/norm(dx);
        adhF1(:,p) = adhF1(:,p) - force;
        adhF2(:,q) = adhF2(:,q) + force;
    end
end