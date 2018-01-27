function flag = IsObject(x,y,r,shape)
% judge if position (x,y) is inside the object
%% for sphere
if strcmp(shape,'Sphere');
    flag = x^2+y^2 < r^2;
end
%% for cube
if strcmp(shape,'Cube');
    flag = 1;
end
%% for rectangle
if strcmp(shape,'Rect');
    flag = y>0;
end
%% for bunny
if strcmp(shape,'Bunny');
    pic = im2bw(imread('shapes\\bunny.png'));
    [a,b] = size(pic);
    flag = ~pic(ceil(a*(r-y)/(2*r)),ceil(b*(x+r)/(2*r)));
end
%% for leaf
if strcmp(shape,'Leaf');
    pic = im2bw(imread('shapes\\leaf.png'));
    [a,b] = size(pic);
    flag = ~pic(ceil(a*(r-y)/(2*r)),ceil(b*(x+r)/(2*r)));
end