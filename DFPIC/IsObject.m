function flag = IsObject(x,y,r,shape)
% judge if position (x,y) is inside the object
switch shape
    case 'Sphere'
        flag = x^2+y^2 < r^2;
    case 'Cube'
        flag = 1;
    case 'Rect'
        flag = y<-1/3*r;
    case 'Dam'
        flag = x<1/4*r;
    case 'Bunny'
        pic = im2bw(imread('shapes\\bunny.png'));
        [a,b] = size(pic);
        flag = ~pic(ceil(a*(r-y)/(2*r)),ceil(b*(x+r)/(2*r)));
    case 'Leaf'
        pic = im2bw(imread('shapes\\leaf.png'));
        [a,b] = size(pic);
        flag = ~pic(ceil(a*(r-y)/(2*r)),ceil(b*(x+r)/(2*r)));
end