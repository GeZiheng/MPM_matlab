function [mi,vi] = ReducedBasis(xi,vi,mp,xp,vp,mode,kernel,dt,iter_time)
for n=1:iter_time
    [~,vp,~,~,Cp] = G2P(xi,vi,xp,vp,vi,1,dt,mode,kernel);
    [mi,vi] = P2G(xp,mp,vp,Cp,xi,mode,kernel);
end