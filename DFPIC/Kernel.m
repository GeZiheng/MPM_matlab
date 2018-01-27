function [base_node,w,dw] = Kernel(x,deg)
% compute 1D B spline weights (linear,quadratic,cubic)
% x is assumed to be scaled in the index space (i.e., it is in a dx=1 grid)
%% Linear Bspline kernel
if deg == 1
    base_node = floor(x);
    w = [ base_node+1-x, x-base_node ];
    dw = [ 1, -1 ];
end

%% Quadratic Bspline kernel
if deg == 2
    base_node = floor(x-0.5);
    w = [ 1/2 * (3/2-x+base_node)^2, 3/4 - (base_node+1-x)^2, 1/2 * (x-base_node-1/2)^2 ];
    dw = [ x-base_node-3/2, 2*(base_node+1-x), x-base_node-1/2 ];
end

%% Cubic Bspline kernel
if deg == 3
    base_node = floor(x)-1;
    w = [ 1/6 * (2-x+base_node)^3, 1/2 * (x-base_node-1)^3 - (x-base_node-1)^2 + 2/3, 1/2 * (base_node+2-x)^3 - (base_node+2-x)^2 + 2/3, 1/6 * (x-base_node-1)^3 ];
    dw = [ -1/2 * (2-x+base_node)^2, 3/2 * (x-base_node-1)^2 - 2 * (x-base_node-1), -3/2 * (base_node+2-x)^2 + 2 * (base_node+2-x), 1/2 * (x-base_node-1)^2 ];
end