function [Fe,Fp,J] = EvolveF(nPts,grad_vp,div_vp,Fe,Fp,J,tc,ts,dt,texture,lambda,mu)
% update deformation gradient
tot_n = 0;
for n = 1:length(nPts)
    for p = tot_n+1:tot_n+nPts(n)
        % evolve Fe
        Ftmp = (eye(2)+dt*grad_vp(:,:,p)) * Fe(:,:,p);
        switch texture(n,:)
            case 'Jello'
                Fe(:,:,p) = Ftmp;
            case 'Snow'
                % check if singular value exceed threshold
                [U,S,V] = svd(Ftmp);
                s = [S(1,1),S(2,2)];
                for p=1:2
                    if s(p) < 1-tc
                        s(p) = 1-tc;
                    end
                    if s(p) > 1+ts
                        s(p) = 1+ts;
                    end
                end
                R = diag(s);
                T = diag(1./s);
                % adjust plasticity
                Fe(:,:,p) = U * R * V';
                Fp(:,:,p) = V * T * S * V' * Fp(:,:,p);
            case 'Water'
                J(p) = (1 + dt*div_vp(p)) * J(p);
            case 'Sand'
                [U,S,V] = svd(Ftmp);
                s = [S(1,1),S(2,2)];
                eps = log(s);
                eps_new = eps - sum(eps)/2;
                delta = norm(eps) + (1+lambda(n)/mu(n)) * sum(eps) * 0.5;
                if delta <= 0   % case 1
                    s_new = s;
                else if sum(eps_new)>0    % case 2
                        s_new = [1,1];
                    else        % case 3
                        Hp = eps - delta * eps_new / norm(eps_new);
                        s_new = exp(Hp);
                    end
                end
                Fe(:,:,p) = U * diag(s_new) * V';
                Fp(:,:,p) = V * diag(s_new./s) * V' * Fp(:,:,p);
        end
    end
    tot_n = p;
end