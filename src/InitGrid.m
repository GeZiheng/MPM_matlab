function [xi,vi,di] = InitGrid(cube,h)
% initialize grid
[gx,gy] = meshgrid(cube(1):h:cube(2),cube(3):h:cube(4));
xi(1,:,:) = gx';
xi(2,:,:) = gy';
vi = zeros(size(xi));
[~,m,n] = size(xi);
di = zeros(m,n);
for i = 1:m
	for j = 1:n
		% compute distance field
		di(i,j) = min([xi(1,i,j)-cube(1),cube(2)-xi(1,i,j),xi(2,i,j)-cube(3),cube(4)-xi(2,i,j)]);
	end
end