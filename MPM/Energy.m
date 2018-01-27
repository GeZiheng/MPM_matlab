function [Ek,Ee,Eg] = Energy(nPts,mp,xp,vp,div_vp,Vp,Fe,J,texture,mu,lambda,g)
% compute energy
Ek = 0;
Ee = 0;
Eg = 0;
tot_n = 0;
for n = 1:length(nPts)
    for k = tot_n+1:tot_n+nPts(n)
        % kinetic energy
        Ek = Ek + 1/2 * mp(k) * vp(:,k)' * vp(:,k);
        % elastic energy
        switch texture(n,:)
            case 'Jello'
                F = Fe(:,:,k);
                [U,S,V] = svd(F);
                R = U*V';
                Psi = mu(n) * norm(F-R)^2 + lambda(n)/2 * (det(F)-1)^2;
            case 'Water'
                Psi = div_vp(k)^2;
            case 'Snow'
                F = Fe(:,:,k);
                [U,S,V] = svd(F);
                R = U*V';
                Psi = mu(n) * norm(F-R)^2 + lambda(n)/2 * (det(F)-1)^2;
            case 'Sand'
                F = Fe(:,:,k);
                [U,S,V] = svd(F);
                R = U*V';
                Psi = mu(n) * norm(F-R)^2 + lambda(n)/2 * (det(F)-1)^2;
        end
        Ee = Ee + Psi * Vp(:,k);
        % gravity energy
        Eg = Eg + mp(k) * xp(2,k) * g;
    end
    tot_n = k;
end