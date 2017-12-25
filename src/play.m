% clear all;
%% parameters
kernel = 2;                      % select kernel: linear, quadratic, or cubic
transfer = 'DFPIC';               % transfer scheme: ppic or dfpic
index = 1:7;                     % index for velocity basis
mu = 5e4;                        % bulk modulus
lambda = 7;                      % stretch coeff
rho = 1000;                      % density
g = 10;                          % gravity acceleration
dt = 1e-4;                       % length of each time step
frame_dt = 0.01;                 % time gap between frames
dx = 0.02;                       % distance between grid nodes
N_step = 10000;                  % number of time steps
iter = 0;                        % number of extra iterations
%% initialization
cube = [0,1,0,1];
[xi,vi,di] = InitGrid(cube,dx);
% water ball dropping
%[xp,vp,mp,nPts] = AddPts([0.5;0.5],0.1,[0;-2],0,rho,'Cube',dx);
% dam break
[xp,vp,mp,nPts] = AddPts([0.3;0.3],0.26,[0;0],0,rho,'Dam',dx);
Sp = zeros(length(index),nPts);            % matrices Sp
Vp = zeros(1,nPts);                        % volume Vp
J = ones(1,nPts);                          % Jacobi J
div_vp = zeros(1,nPts);                    % velocity divergence
data = 0;
%% main loop
frame = 0;
for step = 1:N_step
    fprintf('==================== Step %d ================= \n',step);
    [P,L] = Momentum(mp,xp,vp,[0.5;0.5],Sp,dx,2,3);
    fprintf('Linear momentem: [%3d, %3d]\n',P);
    fprintf('Angular momentem: %3d\n',L);
    div = norm(div_vp);
    fprintf('Divergence: %3d\n',div);
    data(step) = div;
    % transfer particles to grid
    [mi,vi] = P2G(xp,mp,vp,Sp,xi,kernel,transfer,index);
    % compute particle volumes
    if step==1
        Vp = Volume(mi,xi,mp,xp,kernel);
    end
    % compute grid forces
    fi = GridForces(xi,xp,J,Vp,mu,lambda,kernel);
    % update grid velocity
    vi = GridVelocity(mi,vi,fi,di,g,dx,dt);
    % reduced basis iteration
%     for sub_step = 1:iter
%     	[vp,~,Sp] = G2P(xi,vi,xp,kernel,transfer,index);
%     	[~,vi] = P2G(xp,mp,vp,Sp,xi,kernel,transfer,index);
%     end
    % update deformation divergence J
    J = J + div_vp.*J*dt;
    % transfer grid back to particles
    [vp,div_vp,Sp] = G2P(xi,vi,xp,kernel,transfer,index);
    % particle advection
    xp = Advection(xp,vp,xi,dt);
    % draw particles and grid on screen
    if mod( step, floor(frame_dt/dt) ) == 1
        frame = frame + 1;
        rdata = min( max(-div_vp,0) / 20, 1 );
        gdata = min( max(div_vp,0) / 20, 1 );
        bdata = 1 - rdata - gdata;
        sca = scatter(xp(1,:),xp(2,:),5.);
        set(sca,'CData',[rdata;gdata;bdata]');
        hold off;
        grid on;
        axis(cube);
        axis square;
        set(gca,'xtick',cube(1):dx:cube(2),'ytick',cube(3):dx:cube(4),'GridLineStyle','-');
        set(gca,'xticklabel',[],'yticklabel',[]);
        set(gcf,'Position',[300,100,600,600]);
        drawnow
        saveas(figure(1),strcat('./tmp/frame',num2str(frame,'%03d'),'.png'));
    end
end
figure;
plot(data);
