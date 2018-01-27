function v = GetVFromS(vp,Sp,xi,xc,xp,h,transfer,index)
% update coefficient matrix Sp for polypic
r = xi - xp;
[gx,dgx] = CalcFuncG(r(1),xp(1)-xc(1),h);
[gy,dgy] = CalcFuncG(r(2),xp(2)-xc(2),h);
switch transfer
    case 'PPIC'
        basis = [ r(1)*eye(2), r(2)*eye(2), r(1)*r(2)*eye(2), gx*eye(2), gy*eye(2), gx*r(2)*eye(2), r(1)*gy*eye(2), gx*gy*eye(2) ];
    case 'DFPIC'
        [hx,dhx] = CalcFuncH(r(1),xp(1)-xc(1),h);
        [hy,dhy] = CalcFuncH(r(2),xp(2)-xc(2),h);
        basis = [ r(2), 0, gy, 0, r(1), hx, -r(1)*dhy; 
            0, r(1), 0, gx, -r(2), -dhx*r(2), hy; ];
        %basis = [ r(2), 0, gy, 0, gx, -r(1)*dgy; 
            %0, r(1), 0, gx, -dgx*r(2), gy; ];
end
v = vp + basis(:,index) * Sp;