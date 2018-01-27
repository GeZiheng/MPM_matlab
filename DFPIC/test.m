clear all;
%% parameters
kernel = 2;                      % select kernel: linear, quadratic, or cubic
transfer = 'DFPIC';              % transfer scheme: ppic or dfpic
index = 1:7;                     % index for velocity basis
mu = 1e4;                        % bulk modulus
lambda = 7;                      % stretch coeff
rho = 1000;                      % density
g = 0;                           % gravity acceleration
dt = 1e-4;                       % length of each time step
frame_dt = 0.005;                % time gap between frame
dx = 0.1;                        % distance between grid nodes
N_step = 50;                     % number of time steps
%% initialization
cube = [0,1,0,1];
[xi,vi,di] = InitGrid(cube,dx);
[xp,vp,mp,nPts] = quad_velo([0.5;0.5],0.2,dx,rho,'div',100,10);
Sp = zeros(length(index),nPts);   % matrices Sp
data = 0;
%% main loop
for step = 1:N_step
    fprintf('==================== Step %d ================= \n',step);
    % transfer particles to grid
    [mi,vi] = P2G(xp,mp,vp,Sp,xi,kernel,transfer,index);
    quiver(xi(1,:,:),xi(2,:,:),vi(1,:,:),vi(2,:,:),'LineWidth',2,'color','green');
    grid on;
    axis(cube);
    axis square;
    set(gca,'xtick',cube(1):dx:cube(2),'ytick',cube(3):dx:cube(4),'GridLineStyle','-');
    [vp,div_vp,Sp] = G2P(xi,vi,xp,kernel,transfer,index);
    div = norm(div_vp);
    data(step) = div;
    fprintf('Divergence: %3d\n',div);
    drawnow;
    saveas(figure(1),strcat('./tmp/frame',num2str(step,'%03d'),'.png'));
end
figure;
plot(data);