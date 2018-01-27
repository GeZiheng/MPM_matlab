function [xp,vp,mp,nPts,Cp,Vp,Fe,Fp,J,grad_vp,div_vp] = AddSpring(xp,vp,mp,nPts,Cp,Vp,Fe,Fp,J,grad_vp,div_vp,c,r,v,rho,h)
% increase particles to a spring
sample = poissonDisc([512,512],h*60/r);      % sample particles
sample = sample * 2*r/512 - r;
old_nPts = nPts;
for p = 1:length(sample(:,1))
    x = sample(p,1);
    y = sample(p,2);
    if y > 0 || y < -h
		continue;
	end
    nPts = nPts + 1;
	mp(nPts) = rho * (h^2/10);
    xp(1,nPts) = c(1) + x;
    xp(2,nPts) = c(2) + y;
    vp(1,nPts) = v(1);
    vp(2,nPts) = v(2);
end
new_nPts = nPts - old_nPts;
Cp = cat(3,Cp,zeros(2,2,new_nPts));
Vp = [Vp,zeros(1,new_nPts)];
J = [J,ones(1,new_nPts)];
grad_vp = cat(3,grad_vp,zeros(2,2,new_nPts));
div_vp = [div_vp,zeros(1,new_nPts)];
for k=old_nPts+1:nPts
    Fe(:,:,k) = eye(2);     % matrices for elasticity
    Fp(:,:,k) = zeros(2);   % matrices for plasticity
end