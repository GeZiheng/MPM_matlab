clear all;
%% parameters
MODE = 'APIC';                  % select mode: APIC or FLIP
KERNEL = 2;                     % select kernel: linear, quadratic, or cubic
SCHEME = 'Explicit';            % select scheme: explicit or implicit
MODEL = 'None';                 % select model for surface tension: CSS, CSF, or SPH
SAMPLE = 'Uniform';             % select sample mode: uniform or poisson disk
flip = 1;                       % parameter for FLIP & PIC blending
g = 10;                         % gravity acceleration
dt = 5e-5;                      % length of each time step
frame_dt = 0.01;                % time gap between frame
dx = 1/64;                      % distance between grid nodes
N_step = 10000;                 % number of time steps
beta = 10;                      % parameter for adhesion
width = 0.1;
spring_rate = 0.3;
spring_dx = dx;
%% initialization
boundary_cube = [0,1,0,1];
collide_cube = [0.5,0.5,0.08];
[xi,vi,di] = InitGrid(boundary_cube,collide_cube,dx);
nPts = 0;
% water dropping on a leaf
%[xp1,vp1,mp1,nPts(1)] = AddPts([0.35,0.8],0.075,[0,0],0,1000,'Cube',dx);
%[xp2,vp2,mp2,nPts(2)] = AddPts([0.3,0.4],0.3,[0,0],0,400,'Leaf',dx);
% water balls collision
%[xp1,vp1,mp1,nPts1] = AddPts([0.5,0.35],0.1,[0,2],0,1000,'Sphere',dx);
%[xp2,vp2,mp2,nPts2] = AddPts([0.5,0.65],0.1,[0,-2],0,1000,'Sphere',dx);
% water standing on jello
%[xp1,vp1,mp1,nPts(1)] = AddPts([0.5,0.6],0.2,[0,0],0,1000,'Sphere',dx);
%[xp2,vp2,mp2,nPts(2)] = AddPts([0.5,0.45],0.4,[0,0],0,400,'Rect',dx);
% water dripping
%[xp1,vp1,mp1,nPts] = AddPts([0.5,1-width],width,[0,-spring_rate],0,1000,'Rect',dx);
% dam break
[xp1,vp1,mp1,nPts] = AddPts([0.5,0.8],0.1,[0,0],0,400,'Rect',dx,SAMPLE);
texture = 'Sand';
%texture = ['Water';'Jello'];
[mu,lambda,rho,tc,ts,ksi,gamma] = SetPara(texture);
xp = xp1;
vp = vp1;
mp = mp1;
%xp = [xp1,xp2];
%vp = [vp1,vp2];
%mp = [mp1,mp2];
%nPts = nPts1+nPts2;
Vp = zeros(1,sum(nPts));     % volume Vp
J = ones(1,sum(nPts));       % Jacobi J
grad_vp = zeros(2,2,sum(nPts)); % velocity gradient
div_vp = zeros(1,sum(nPts));    % velocity divergence
data = 0;
Fe = [];
Fp = [];
for k=1:sum(nPts)
    Fe(:,:,k) = eye(2);     % matrices for elasticity
    Fp(:,:,k) = zeros(2);   % matrices for plasticity
end
set(gcf,'Position',get(0,'ScreenSize'))
%% main loop
frame = 0;
for step = 1:N_step
    % water dripping
%      if max(xp(2,:))<=1-0.8*spring_dx
%          [xp,vp,mp,nPts,Cp,Vp,Fe,Fp,J,grad_vp,div_vp] = AddSpring(xp,vp,mp,nPts,Cp,Vp,Fe,Fp,J,grad_vp,div_vp,[0.5,1.0],width,[0,-spring_rate],1000,spring_dx);
%      end
    fprintf('==================== Step %d ================= \n',step);
    % transfer particles to grid
    [mi,vin] = P2G(xp,mp,vp,grad_vp,xi,MODE,KERNEL);
    % compute particle volumes
    if step==1
        Vp = Volume(mi,xi,mp,xp,KERNEL);
    end
    % compute grid forces
    fi = GridForces(nPts,xi,xp,mp,Fe,Fp,J,Vp,mu,lambda,ksi,gamma,KERNEL,texture,MODEL,beta);
    % update grid velocity
    vi = GridVelocity(mi,vin,fi,di,g,dx,dt,SCHEME);
    % reduced basis iteration
    %[mi,vi] = ReducedBasis(xi,vi,mp,xp,vp,MODE,KERNEL,dt,10);
    % update deformation gradient
    [Fe,Fp,J] = EvolveF(nPts,grad_vp,div_vp,Fe,Fp,J,tc,ts,dt,texture,lambda,mu);
    % transfer grid back to particles
    [xp,vp,grad_vp,div_vp] = G2P(xi,vi,xp,vp,vin,flip,dt,MODE,KERNEL);
    % draw particles and grid on screen
    if mod( step, floor(frame_dt/dt) ) == 1
        frame = frame + 1;
        DrawPts(xp,nPts,texture);
        fill([collide_cube(1)-collide_cube(3),collide_cube(1),collide_cube(1)+collide_cube(3),collide_cube(1)],[collide_cube(2),collide_cube(2)-collide_cube(3),collide_cube(2),collide_cube(2)+collide_cube(3)],'b')
        hold off;
        grid on;
        axis(boundary_cube);
        axis square;
        set(gca,'xtick',boundary_cube(1):dx:boundary_cube(2),'ytick',boundary_cube(3):dx:boundary_cube(4),'GridLineStyle','-');
        set(gca,'xticklabel',[],'yticklabel',[]);
        set(gcf,'Position',[300,100,600,600]);
        drawnow
        saveas(figure(1),strcat('./tmp/frame',num2str(frame,'%03d'),'.png'));
    end
    % compute linear & angular momentum
    Dp = dx^2/4;
    [P,L] = Momentum(mp,xp,vp,grad_vp,Dp);
    [Ek,Ee,Eg] = Energy(nPts,mp,xp,vp,div_vp,Vp,Fe,J,texture,mu,lambda,g);
    data(step) = Ee;
    fprintf('Linear momentem: [%3d, %3d]\n',P);
    fprintf('Angular momentem: %3d\n',L);
    fprintf('Energy: [%3d,%3d,%3d]\n',Ek,Ee,Eg);
    fprintf('Total: %3d\n',Ek+Ee+Eg);
end