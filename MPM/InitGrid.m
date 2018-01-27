function [xi,vi,di] = InitGrid(boundary_cube,collide_cube,h)
% initialize grid
[gx,gy] = meshgrid(boundary_cube(1):h:boundary_cube(2),boundary_cube(3):h:boundary_cube(4));
xi(1,:,:) = gx';
xi(2,:,:) = gy';
vi = zeros(size(xi));
[ch,m,n] = size(xi);
di = zeros(m,n);
for i = 1:m
	for j = 1:n
		% compute distance field
		di(i,j) = min([xi(1,i,j)-boundary_cube(1),boundary_cube(2)-xi(1,i,j),xi(2,i,j)-boundary_cube(3),boundary_cube(4)-xi(2,i,j)]);
        di(i,j) = min([di(i,j),abs(xi(1,i,j)-collide_cube(1))+abs(xi(2,i,j)-collide_cube(2))-collide_cube(3)]);
	end
end