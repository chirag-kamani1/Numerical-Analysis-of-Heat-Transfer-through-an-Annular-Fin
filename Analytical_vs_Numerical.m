clear;clc;
k = input('ENTER THE VALUE OF THERMAL CONDUCTIVITY in W/m-°C: ');
h = input('ENTER THE VALUE OF HEAT TRANSFER COFFICIENT in W/m^2-°C : ');
t = 0.02; % Thickness of Annular fin
r1 = 0.05; % Inner Radius of Annular fin
r2 = 0.15; % Outer Radius of Annular fin
dr = 0.001;
r = r1:dr:r2;
nr = length(r);
Tamb = 25;
T_base = input('ENTER THE VALUE OF BASE TEMPERATURE FOR THE FIN IN °C: ');
qo = T_base-Tamb;
% ANALYTICAL SOLUTION
b = sqrt(2*h/(k*t));
Io_i = besseli(0,b*r1);
Ko_i = besselk(0,b*r1);
Io_o = besseli(0,b*r2);
Ko_o = besselk(0,b*r2);
I1_i = besseli(1,b*r1);
K1_i = besselk(1,b*r1);
I1_o = besseli(1,b*r2);
K1_o = besselk(1,b*r2);
c2 = (-qo*(h*Io_o+k*b*I1_o))/(-Ko_i*(h*Io_o+k*b*I1_o)+(h*Ko_o-k*b*K1_o)*Io_i);
c1 = (qo-c2*Ko_i)/Io_i;
Q_analytical = c1*besseli(0,b*r) + c2*besselk(0,b*r);
heat = -k*2*pi*r1*t*(c1*b*besseli(1,b*r1)-c2*b*besselk(1,b*r1));
% NUMERICAL SOLUTION USING FINITITE DIFF. METHOD & GAUSS SEIDEL METHOD
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
tol = 10E-6;
[x] = GaussSiedelMethod(A,B,x0,tol);
Q_numerical = zeros(1,nr);
Q_numerical(1) = qo;
Q_numerical(2:end) = x;
plot(r,Q_analytical,'m','LineWidth',1.5);
hold on;
plot(r,Q_numerical,'o', 'MarkerFaceColor', 'g', 'MarkerSize',4);
xlabel('RADIAL DISTANCE (ri to ro) in m', FontWeight='bold');
ylabel('θ = (T-Tamb) in °C', 'FontWeight','bold');
legend('ANALYTICAL','NUMERICAL');
title('TEMPERATURE DISTRIBUTION FOR ANNULAR FIN ALONG RADIAL DIRECTION');
grid on;
function [x] = GaussSiedelMethod(A,B,x0,tol)
n = length(A);
x = x0;
iterations = 0;
err = max(A);
if(CheckDiagDom(A))
while(err>tol)
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