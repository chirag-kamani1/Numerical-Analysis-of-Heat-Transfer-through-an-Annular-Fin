clear; clc;
t = 0.02;
r1 = 0.05;
r2 = 0.15;
q = 3;
disp("PRESS 1 FOR CASE 1 : EFFECT OF BASE TEMPERATURE");
disp("PRESS 2 FOR CASE 2 : EFFECT OF HEAT TRANSFER COEFFICIENT");
disp("PRESS 3 FOR CASE 3 : EFFECT OF FIN MATERIAL");
disp("PRESS 4 FOR CASE 4 : GRID INDEPENDENCE TEST");
P = zeros(3,4);
choice = input("ENTER YOUR CHOICE: ");
switch (choice)
case 1
n = input("ENTER THE NUMBER OF NODES FOR DISCRETIZATION: ");
h = input("ENTER THE VALUE OF CONVECTIVE HEAT TRANSFER COEFFICIENT: ");
k = input("ENTER THE VALUE OF THERMAL CONDUCTIVITY OF FIN MATERIAL: ");
Tamb = input("ENTER THE AMBIENT TEMPERATURE IN CELSIUS: ");
dr = (r2-r1)/(n-1);
r = r1:dr:r2;
nr = length(r);
for z = 1:q
T_base = input("ENTER THE BASE TEMPERATURE IN CELSIUS: ");
qo = T_base-Tamb;
P(1,z) = qo;
l = (1/dr^2-0.5./(r.*dr));
m = (-2/dr^2-2*h/(k*t));
n = (1/dr^2+0.5./(r.*dr));
b = sqrt(2*h/(k*t));
A = zeros(nr-1,nr-1);
A(1,1) = m;
A(1,2) = n(2);
for i = 2:nr-2
A(i,i) = m;
A(i,i-1) = l(i+1);
A(i,i+1) = n(i+1);
end
A(nr-1,nr-2) = l(end)+n(end);
A(nr-1,nr-1) = m-n(end)*h*2*dr/k;
B = zeros(nr-1,1);
B(1) = -qo*l(2);
x0 = zeros(nr-1,1);
acc = 10E-6;
[x,iterations] = GaussSiedelMethod(A,B,x0,acc);
Q_numerical = zeros(nr,1);
Q_numerical(1) = qo;
Q_numerical(2:end) = x;
plot(r,Q_numerical,'LineWidth',1.3);
hold on;
xlabel('RADIAL DIRECTION (r1 to r2) in m','FontWeight','bold');
ylabel('θ = (T-Tamb) in °C', 'FontWeight','bold');
title('TEMPERATURE DISTRIBUTION ALONG RADIAL DIRECTION');
grid on;
end
legend('Tbase1','Tbase2','Tbase3');
case 2
n = input("ENTER THE NUMBER OF NODES FOR DISCRETIZATION: ");
k = input("ENTER THE VALUE OF THERMAL CONDUCTIVITY OF FIN MATERIAL: ");
Tamb = input("ENTER THE AMBIENT TEMPERATURE IN CELSIUS: ");
T_base = input("ENTER THE BASE TEMPERATURE IN CELSIUS: ");
dr = (r2-r1)/(n-1);
r = r1:dr:r2;
nr = length(r);
qo = T_base-Tamb;
for z = 1:q
h = input("ENTER THE VALUE OF CONVECTIVE HEAT TRANSFER COEFFICIENT: ");
P(2,z) = h;
l = (1/dr^2-0.5./(r.*dr));
m = (-2/dr^2-2*h/(k*t));
n = (1/dr^2+0.5./(r.*dr));
b = sqrt(2*h/(k*t));
A = zeros(nr-1,nr-1);
A(1,1) = m;
A(1,2) = n(2);
for i = 2:nr-2
A(i,i) = m;
A(i,i-1) = l(i+1);
A(i,i+1) = n(i+1);
end
A(nr-1,nr-2) = l(end)+n(end);
A(nr-1,nr-1) = m-n(end)*h*2*dr/k;
B = zeros(nr-1,1);
B(1) = -qo*l(2);
x0 = zeros(nr-1,1);
acc = 10E-6;
[x,iterations] = GaussSiedelMethod(A,B,x0,acc);
Q_numerical = zeros(nr,1);
Q_numerical(1) = qo;
Q_numerical(2:end) = x;
plot(r,Q_numerical, 'LineWidth',1.3);
hold on;
xlabel('RADIAL DIRECTION (r1 to r2) in m','FontWeight','bold');
ylabel('θ = (T-Tamb) in °C', 'FontWeight','bold');
title('TEMPERATURE DISTRIBUTION ALONG RADIAL DIRECTION');
grid on;
end
legend('h1','h2','h3');
case 3
n = input("ENTER THE NUMBER OF NODES FOR DISCRETIZATION: ");
h= input("ENTER THE VALUE OF CONVECTIVE HEAT TRANSFER COEFFICIENT: ");
Tamb = input("ENTER THE AMBIENT TEMPERATURE IN CELSIUS: ");
T_base = input("ENTER THE BASE TEMPERATURE IN CELSIUS: ");
dr = (r2-r1)/(n-1);
r = r1:dr:r2;
nr = length(r);
qo = T_base-Tamb;
for z = 1:q
k = input("ENTER THE VALUE OF THERMAL CONDUCTIVITY OF FIN MATERIAL: ");
P(3,z) = k;
l = (1/dr^2-0.5./(r.*dr));
m = (-2/dr^2-2*h/(k*t));
n = (1/dr^2+0.5./(r.*dr));
b = sqrt(2*h/(k*t));
A = zeros(nr-1,nr-1);
A(1,1) = m;
A(1,2) = n(2);
for i = 2:nr-2
A(i,i) = m;
A(i,i-1) = l(i+1);
A(i,i+1) = n(i+1);
end
A(nr-1,nr-2) = l(end)+n(end);
A(nr-1,nr-1) = m-n(end)*h*2*dr/k;
B = zeros(nr-1,1);
B(1) = -qo*l(2);
x0 = zeros(nr-1,1);
acc = 10E-6;
[x,iterations] = GaussSiedelMethod(A,B,x0,acc);
Q_numerical = zeros(nr,1);
Q_numerical(1) = qo;
Q_numerical(2:end) = x;
plot(r,Q_numerical, 'LineWidth',1.3);
hold on;
xlabel('RADIAL DIRECTION (r1 to r2) in m','FontWeight','bold');
ylabel('θ = (T-Tamb) in °C', 'FontWeight','bold');
title('TEMPERATURE DISTRIBUTION ALONG RADIAL DIRECTION');
grid on;
end
legend('K1','K2','K3');
case 4
k = input("ENTER THE VALUE OF THERMAL CONDUCTIVITY OF FIN MATERIAL: ");
h= input("ENTER THE VALUE OF CONVECTIVE HEAT TRANSFER COEFFICIENT: ");
Tamb = input("ENTER THE AMBIENT TEMPERATURE IN CELSIUS: ");
T_base = input("ENTER THE BASE TEMPERATURE IN CELSIUS: ");
qo = T_base-Tamb;
for z = 1:q
n = input("ENTER THE NUMBER OF NODES FOR DISCRETIZATION: ");
P(4,z) = n;
dr = (r2-r1)/(n-1);
r = r1:dr:r2;
nr = length(r);
l = (1/dr^2-0.5./(r.*dr));
m = (-2/dr^2-2*h/(k*t));
n = (1/dr^2+0.5./(r.*dr));
b = sqrt(2*h/(k*t));
A = zeros(nr-1,nr-1);
A(1,1) = m;
A(1,2) = n(2);
for i = 2:nr-2
A(i,i) = m;
A(i,i-1) = l(i+1);
A(i,i+1) = n(i+1);
end
A(nr-1,nr-2) = l(end)+n(end);
A(nr-1,nr-1) = m-n(end)*h*2*dr/k;
B = zeros(nr-1,1);
B(1) = -qo*l(2);
x0 = zeros(nr-1,1);
acc = 10E-6;
[x,iterations] = GaussSiedelMethod(A,B,x0,acc);
Q_numerical = zeros(nr,1);
Q_numerical(1) = qo;
Q_numerical(2:end) = x;
plot(r,Q_numerical, 'LineWidth',1.3);
hold on;
xlabel('RADIAL DIRECTION (r1 to r2) in m','FontWeight','bold');
ylabel('θ = (T-Tamb) in °C', 'FontWeight','bold');
title('TEMPERATURE DISTRIBUTION ALONG RADIAL DIRECTION');
grid on;
end
legend('n1','n2','n3');
otherwise
fprintf('Invalid Input..Try Again.\n');
return;
end
function [x,iterations]=GaussSiedelMethod(A,B,x0,acc)
n = length(A);
x = x0;
iterations = 0;
err = max(A);
if(CheckDiagDom(A))
while(err>acc)
xdummy = x;
iterations = iterations+1;
for i = 1:n
sumOfOther = sum(A(i,1:n)*x(1:n))-A(i,i)*x(i);
x(i) = (B(i)-sumOfOther)/A(i,i);
end
err = sqrt(sum((x-xdummy).^2)/n);
end
else
iterations = 0;
end
end
function isDiagDom = CheckDiagDom(M)
isDiagDom = true;
for i = 1:size(M)
sumofOffDiag = sum(abs(M(1,:)))-abs(M(i,i));
if (abs(M(i,i))<sumofOffDiag)
isDiagDom = false;
break;
end
end
end