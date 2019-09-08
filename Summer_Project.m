%% 玡矗

clear

%% Parameters

% m-kg-N-s-Pa-rad
% r1, r2ゼ
L1 = 9.14; % m
L2 = L1/cos(pi/4.0); % m
E = 200.0e9; % Pa
rho = 7860.0; % kg/m^3
sigy = 250.0e6; % Pa


%% 程ㄎて

% fmincon setup
x0 = [0.1;0.05];
A = [];
b = [];
Aeq = [];
beq = [];
vlb = [0.001;0.001];
vub = [0.5;0.5];
options = optimoptions('fmincon');

[x,fval,exitflag] = fmincon(@mincon,x0,A,b,Aeq,beq,vlb,vub,@limit,options);

%% fmincon Functions

function f = mincon(r)
L1 = 9.14;
L2 = L1/cos(pi/4.0);
A1 = r(1)*r(1)*pi;
A2 = r(2)*r(2)*pi;
rho = 7860.0;

% m1 + m2
f = A1*L1*rho*6.0 + A2*L2*rho*4.0;

end

function [g,ceq] = limit(r)
% K & Q calculation
% 爹秆DOF

L1 = 9.14;
L2 = L1/cos(pi/4.0);
A1 = r(1)*r(1)*pi;
A2 = r(2)*r(2)*pi;
E = 200.0e9;
sigy = 250.0e6;

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

% KQ = F
% Q(1,7:10) = 0
Q = kmat(1:8,1:8)\[0;0;0;-10000000;0;0;0;-10000000];
Q = [Q;0;0;0;0];

%―sigma
sig1 = sig(pi,E,L1,Q,5,6,9,10);
sig2 = sig(pi,E,L1,Q,1,2,5,6);
sig3 = sig(pi,E,L1,Q,7,8,11,12);
sig4 = sig(pi,E,L1,Q,3,4,7,8);
sig5 = sig(pi*3/2,E,L1,Q,5,6,7,8);
sig6 = sig(pi*3/2,E,L1,Q,1,2,3,4);
sig7 = sig(pi*3/4,E,L2,Q,7,8,9,10);
sig8 = sig(pi*5/4,E,L2,Q,5,6,11,12);
sig9 = sig(pi*3/4,E,L2,Q,3,4,5,6);
sig10 = sig(pi*5/4,E,L2,Q,1,2,7,8);

% limitations matrix
% g = [s1-sigy; s2-sigy; s3-sigy; s4-sigy; s5-sigy; s6-sigy; s7-sigy; s8-sigy; s9-sigy; s10-sigy; sqrt(Q(3)*Q(3)+Q(4)*Q(4))-0.02];
s = [sig1,sig2,sig3,sig4,sig5,sig6,sig7,sig8,sig9,sig10];
max_stress = max(abs(s));
g = [max_stress-sigy;sqrt(Q(3)^2+Q(4)^2)-0.02];
ceq = [];
end

%% Calculations function
function k = ktheta(theta,A,E,L)
% K matrix calculation for each truss
% 计node计タ
k = A*E/L*[cos(theta)^2 cos(theta)*sin(theta) -cos(theta)^2 -cos(theta)*sin(theta);
    cos(theta)*sin(theta) sin(theta)^2 -cos(theta)*sin(theta) -sin(theta)^2;
    -cos(theta)^2 -cos(theta)*sin(theta) cos(theta)^2 cos(theta)*sin(theta);
    -cos(theta)*sin(theta) -sin(theta)^2 cos(theta)*sin(theta) sin(theta)^2];
end

function big = ktrans(M,d1,d2,d3,d4)
%盢4X4痻皚耎12X12痻皚よ獽
big = zeros(12);
% big(d1,d1) = M(1,1);
% big(d1,d2) = M(1,2);
% big(d1,d3) = M(1,3);
% big(d1,d4) = M(1,4);
% 
% big(d2,d1) = M(2,1);
% big(d2,d2) = M(2,2);
% big(d2,d3) = M(2,3);
% big(d2,d4) = M(2,4);
% 
% big(d3,d1) = M(3,1);
% big(d3,d2) = M(3,2);
% big(d3,d3) = M(3,3);
% big(d3,d4) = M(3,4);
% 
% big(d4,d1) = M(4,1);
% big(d4,d2) = M(4,2);
% big(d4,d3) = M(4,3);
% big(d4,d4) = M(4,4);
for i = 1:4
    im = [d1 d2 d3 d4];
    for j = 1:4
        jm = [d1 d2 d3 d4];
        big(im(i),jm(j)) = M(i,j);
    end
end
end

function f = sig(theta,E,L,Q,d1,d2,d3,d4)
% sigma calculation for each truss
% 计node计タ
f = E/L*[-cos(theta) -sin(theta) cos(theta) sin(theta)]...
    *[Q(d1);Q(d2);Q(d3);Q(d4)];
end
