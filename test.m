clear

% m-kg-N-s-Pa-rad
r1 = 0.1; % m
r2 = 0.05; % m
A1 = r1^2*pi;
A2 = r2^2*pi;
L1 = 9.14; % m
L2 = L1/cos(pi/4.0); % m
E = 200.0e9; % Pa
rho = 7860.0; % kg/m^3
sigy = 250.0e6; % Pa

% k matrix calculation
k1 = ktheta(pi,A1,E,L1); % 5 6 9 10
k2 = ktheta(pi,A1,E,L1); % 1 2 5 6
k3 = ktheta(pi,A1,E,L1); % 7 8 11 12
k4 = ktheta(pi,A1,E,L1); % 3 4 7 8
k5 = ktheta(pi*3/2,A1,E,L1); % 5 6 7 8
k6 = ktheta(pi*3/2,A1,E,L1); % 1 2 3 4
k7 = ktheta(pi*3/4,A2,E,L2); % 7 8 9 10
k8 = ktheta(pi*5/4,A2,E,L2); % 5 6 11 12
k9 = ktheta(pi*3/4,A2,E,L2); % 3 4 5 6
k10 = ktheta(pi*5/4,A2,E,L2); % 1 2 7 8

% mix all 4X4 matrix into 12X12
kmat = ktrans(k1,5,6,9,10) + ktrans(k2,1,2,5,6) + ...
    ktrans(k3,7,8,11,12) + ktrans(k4,3,4,7,8) + ...
    ktrans(k5,5,6,7,8) + ktrans(k6,1,2,3,4) + ...
    ktrans(k7,7,8,9,10) + ktrans(k8,5,6,11,12) + ...
    ktrans(k9,3,4,5,6) + ktrans(k10,1,2,7,8);

function k = ktheta(theta,A,E,L)
% K matrix calculation for each truss
% 數字小的node到數字大的為正
k = A*E/L*[cos(theta)^2 cos(theta)*sin(theta) -cos(theta)^2 -cos(theta)*sin(theta);
    cos(theta)*sin(theta) sin(theta)^2 -cos(theta)*sin(theta) -sin(theta)^2;
    -cos(theta)^2 -cos(theta)*sin(theta) cos(theta)^2 cos(theta)*sin(theta);
    -cos(theta)*sin(theta) -sin(theta)^2 cos(theta)*sin(theta) sin(theta)^2];
end

function big = ktrans(M,d1,d2,d3,d4)
%將4X4矩陣擴充為12X12矩陣方便相加
big = zeros(12);
for i = 1:4
    im = [d1 d2 d3 d4];
    for j = 1:4
        jm = [d1 d2 d3 d4];
        big(im(i),jm(j)) = M(i,j);
    end
end
end