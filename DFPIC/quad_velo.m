function [xp,vp,mp,nPts,M,I] = quad_velo(xc,r,h,rho,mode,inte,freq)
% initialize particles
sample = poissonDisc([512,512],h*60/r);      % sample particles
sample = sample * 2*r/512 - r;
[nPts,~] = size(sample);
xp = sample' + xc*ones(1,nPts);
switch mode
    case 'dilation'
        vp = inte*sample';
    case 'collision'
        vp(1,:) = inte*sample(:,1)';
        vp(2,:) = -inte*sample(:,2)';
    case 'shear'
        vp(1,:) = inte*sample(:,2)';
        vp(2,:) = inte*zeros(1,nPts);
    case 'rotation'
        vp(1,:) = -inte*sample(:,2)';
        vp(2,:) = inte*sample(:,1)';
    case 'random'
        vp = inte*(rand(2,nPts)-1/2);
    case 'div'
        % 0 rotation
        vp(1,:) = inte*cos(freq*sample(:,1)').*sin(freq*sample(:,2)');
        vp(2,:) = inte*sin(freq*sample(:,1)').*cos(freq*sample(:,2)');
    case 'rot'
        % 0 divergence
        vp(1,:) = inte*cos(freq*sample(:,1)').*sin(freq*sample(:,2)');
        vp(2,:) = -inte*sin(freq*sample(:,1)').*cos(freq*sample(:,2)');
    case 'divrot'
        % constant divergence and constant rotation
        vp(1,:) = inte*( sin(freq*sample(:,1)').*cos(freq*sample(:,2)') + cos(freq*sample(:,1)').*sin(freq*sample(:,2)') );
        vp(2,:) = inte*( cos(freq*sample(:,1)').*sin(freq*sample(:,2)') - sin(freq*sample(:,1)').*cos(freq*sample(:,2)') );
end
mp = ones(1,nPts) * rho * (h^2/10);  % mass of particles
M = sum(mp);
I = 0;          % calculate rotational inertia
for p = 1:nPts
    I = I + mp(p)*norm(xp(:,p)-xc)^2;
end