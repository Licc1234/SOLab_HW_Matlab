clear

global r1
r1 = 0.1;
global r2
r2 = 0.05;
global L1
L1 = 9.14;
global L2
L2 = L1/cos(pi/4);
global A1;
A1 = r1*r1*pi;
global A2;
A2 = r2*r2*pi;
global E
E = 200*(10^9);
global rho
rho = 7860;
global sigy
sigy = 250*(10^6);
global m1;
m1 = A1*L1*rho*6;
global m2;
m2 = A2*L2*rho*6;

% K calculation
% 註解為DOF
k1 = krad1(pi); % 5 6 9 10
k2 = krad1(pi); % 1 2 5 6
k3 = krad1(pi); % 7 8 11 12
k4 = krad1(pi); % 3 4 7 8
k5 = krad1(pi*3/2); % 5 6 7 8
k6 = krad1(pi*3/2); % 1 2 3 4
k7 = krad2(pi*3/4); % 7 8 9 10
k8 = krad2(pi*5/4); % 5 6 11 12
k9 = krad2(pi*3/4); % 3 4 5 6
k10 = krad2(pi*5/4); % 1 2 7 8
kmat = ktrans(k1,5,6,9,10) + ktrans(k2,1,2,5,6) + ktrans(k3,7,8,11,12) + ktrans(k4,3,4,7,8) + ktrans(k5,5,6,7,8) + ktrans(k6,1,2,3,4) + ktrans(k7,7,8,9,10) + ktrans(k8,5,6,11,12) + ktrans(k9,3,4,5,6) + ktrans(k10,1,2,7,8);
kmat2 = kmat(1:8, 1:8);
% KQ = F
Q = inv(kmat2)*[0;0;0;-10000000;0;0;0;-10000000];
Q = [Q; 0; 0; 0; 0];

function k = krad1(r)
%計算短桿
%數字小的node到數字大的為正
global A1;
global E;
global L1;
k = A1*E/L1*[cos(r).*cos(r) cos(r).*sin(r) -cos(r).*cos(r) -cos(r).*sin(r); cos(r).*sin(r) sin(r).*sin(r) -cos(r).*sin(r) -sin(r).*sin(r); -cos(r).*cos(r) -cos(r).*sin(r) cos(r).*cos(r) cos(r).*sin(r); -cos(r).*sin(r) -sin(r).*sin(r) cos(r).*sin(r) sin(r).*sin(r)];
end

function k = krad2(r)
%計算長桿
%數字小的node到數字大的為正
global A2;
global E;
global L2;
k = A2*E/L2*[cos(r).*cos(r) cos(r).*sin(r) -cos(r).*cos(r) -cos(r).*sin(r); cos(r).*sin(r) sin(r).*sin(r) -cos(r).*sin(r) -sin(r).*sin(r); -cos(r).*cos(r) -cos(r).*sin(r) cos(r).*cos(r) cos(r).*sin(r); -cos(r).*sin(r) -sin(r).*sin(r) cos(r).*sin(r) sin(r).*sin(r)];
end

function big = ktrans(M,d1,d2,d3,d4)
%將4X4矩陣擴充為12X12矩陣方便相加
big = zeros(12);
big(d1,d1) = M(1,1);
big(d1,d2) = M(1,2);
big(d1,d3) = M(1,3);
big(d1,d4) = M(1,4);
big(d2,d1) = M(2,1);
big(d2,d2) = M(2,2);
big(d2,d3) = M(2,3);
big(d2,d4) = M(2,4);
big(d3,d1) = M(3,1);
big(d3,d2) = M(3,2);
big(d3,d3) = M(3,3);
big(d3,d4) = M(3,4);
big(d4,d1) = M(4,1);
big(d4,d2) = M(4,2);
big(d4,d3) = M(4,3);
big(d4,d4) = M(4,4);
end