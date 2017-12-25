function Sp = EvolveS(Sp,wip,vi,xi,xc,xp,h,transfer,index)
% update coefficient matrix Sp for polypic
r = xi - xp;
q = xp - xc;
b = h^2 / 4;
dx = ( h^2 - 4*q(1)^2 )^2 * ( 3*h^2 - 4*q(1)^2 ) / ( 16 * h^2 );
dy = ( h^2 - 4*q(2)^2 )^2 * ( 3*h^2 - 4*q(2)^2 ) / ( 16 * h^2 );
[gx,dgx] = CalcFuncG(r(1),q(1),h);
[gy,dgy] = CalcFuncG(r(2),q(2),h);
switch transfer
    case 'PPIC'
        c = b^2;
        ex = b * dx;
        ey = b * dy;
        f = dx * dy;
        coeff = [ r(1)/b * eye(2);
            r(2)/b * eye(2);
            r(1)*r(2)/c * eye(2); 
            gx/dx * eye(2); 
            gy/dy * eye(2); 
            gx*r(2)/ex * eye(2); 
            r(1)*gy/ey * eye(2); 
            gx*gy/f * eye(2);
            ];
    case 'DFPIC'
        c = 2*b;
        [hx,dhx] = CalcFuncH(r(1),q(1),h);
        [hy,dhy] = CalcFuncH(r(2),q(2),h);
        exx = ( 7*h^6 - 26*h^4*q(1)^2 + 64*h^2*q(1)^4 - 32*q(1)^6 ) / ( 16*h^2 );
        exy = ( h^2 - 4*q(1)^2 ) * q(1) * ( h^2 - 4*q(2)^2 ) * q(2) / ( 8*h^2 );
        eyy = ( 7*h^6 - 26*h^4*q(2)^2 + 64*h^2*q(2)^4 - 32*q(2)^6 ) / ( 16*h^2 );
        %ex = ( 7*h^4 - 24*h^2*q(1)^2 + 48*q(1)^4 ) / 16;
        %ey = ( 7*h^4 - 24*h^2*q(2)^2 + 48*q(2)^4 ) / 16;
        E = [ exx, exy;
             exy, eyy; ];
        F = [ hx, -dhx*r(2);
             -r(1)*dhy, hy; ];
        coeff = [ r(2)/b, 0; 
            0, r(1)/b; 
            gy/dy, 0; 
            0, gx/dx; 
            r(1)/c, -r(2)/c;
            E\F;
            ];
end
Sp = Sp + wip*coeff(index,:)*vi;