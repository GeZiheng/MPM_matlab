function [mu,lambda,rho,tc,ts,ksi,gamma] = SetPara(texture)
% set elasticity parameters
for n = 1:size(texture,1)
    switch texture(n,:)
        case 'Jello'
            E = 1e6;                                % Young's modulus
            nu = 0.3;                                 % Poisson's ratio
            mu(n) = E / (2 * (1 + nu));
            lambda(n) = E * nu / ((1 + nu) * (1 - 2 * nu));
            rho(n) = 400;
            tc(n) = 1;                       % critical compression
            ts(n) = 10;                       % critical stretch
            ksi(n) = 0;                       % hardening coefficient
            gamma(n) = 0;
        case 'Water'
            mu(n) = 1e4;
            lambda(n) = 7;
            rho(n) = 1000;
            tc(n) = 1;                       % critical compression
            ts(n) = 10;                       % critical stretch
            ksi(n) = 0;                      % visid coeff
            gamma(n) = 0.05;                   % surface tension coeff
        case 'Snow'
            E(n) = 1.4e5;                                % Young's modulus
            nu(n) = 0.2;                                 % Poisson's ratio
            mu(n) = E / (2 * (1 + nu));
            lambda(n) = E * nu / ((1 + nu) * (1 - 2 * nu));
            rho(n) = 400;
            tc(n) = 2.5e-2;                       % critical compression
            ts(n) = 7.5e-3;                       % critical stretch
            ksi(n) = 10;                          % hardening coefficient
            gamma(n) = 0;
        case 'Sand'
            E(n) = 1.4e5;                                % Young's modulus
            nu(n) = 0.2;                                 % Poisson's ratio
            mu(n) = E / (2 * (1 + nu));
            lambda(n) = E * nu / ((1 + nu) * (1 - 2 * nu));
            rho(n) = 400;
            tc(n) = 2.5e-2;                       % critical compression
            ts(n) = 7.5e-3;                       % critical stretch
            ksi(n) = 10;                          % hardening coefficient
            gamma(n) = 0;
    end
end