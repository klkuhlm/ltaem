
% http://www.maa.org/external_archive/joma/Volume8/Kalman/General.html

e0 = [6.0,2.0; 
      1.5,4.0];  % (x,y) coordinate shift
alfa = [pi/8; 
	3*pi/4];
ax = [3.0,0.75; 
      2.0,1.0];  % (a,b) major & minor axes

% plotting related (elliptical coordinates)
% --------------------------------------------------

% eccentricity
ec = [sqrt(1 - ax(1,1)^2/ax(1,2)^2);
      sqrt(1 - ax(2,1)^2/ax(2,2)^2)];

% semi-focal radius (f)
f = [ec(1)*ax(1,1);
     ec(2)*ax(2,1)];

% elliptical radii
eta = [acosh(1/ec(1));
       acosh(1/ec(2))]; 

% equispaced angular variable (
phi = linspace(-pi,pi,100)*1i;

z1 = complex(e0(1,1),e0(1,2)) + f(1)*cosh(eta(1) + phi)*exp(alfa(1)*1i);
z2 = complex(e0(2,1),e0(2,2)) + f(2)*cosh(eta(2) + phi)*exp(alfa(2)*1i);

figure(1);
plot(real(z1),imag(z1),'r-');
hold on;
plot(real(z2),imag(z2),'g-');
grid on;
axis equal;

% --------------------------------------------------

% http://www.mathworks.com/matlabcentral/fileexchange/28318-conics-intersection

% The homogeneous representation of a conic is a matrix
% m = [A C D; C B E; D E F] that represents the equation 
%  A x^2 + B y^2 + 2C xy + 2D x + 2Ey + F = 0

% setup canonical matrix representation of ellipse 1
a = ax(1,1); b = ax(1,2);

A = cos(alfa(1))^2/a^2 + sin(alfa(1))^2/b^2;
B = sin(alfa(1))^2/a^2 + cos(alfa(1))^2/b^2;
C = -cos(alfa(1))*sin(alfa(1))*(1/a^2 - 1/b^2);

x0 = e1(1); y0 = e1(2);

D = -(2*x0*A + y0*C)/2;
E = -(2*y0*B + x0*C)/2;
F = A*x0^2 + C*x0*y0 + B*y0^2 - 1;

Q1 = [A, C, D;
      C, B, E;
      D, E, F]

figure(2)
plotConic(Q1,'ro');

% setup canonical matrix representation of ellipse 2
a = ax(2,1); b = ax(2,2);

A = cos(alfa(2))^2/a^2 + sin(alfa(2))^2/b^2;
B = sin(alfa(2))^2/a^2 + cos(alfa(2))^2/b^2;
C = -cos(alfa(2))*sin(alfa(2))*(1/a^2 - 1/b^2);

x0 = e2(1); y0 = e2(2);

D = -(2*x0*A + y0*C)/2;
E = -(2*y0*B + x0*C)/2;
F = A*x0^2 + C*x0*y0 + B*y0^2 - 1;

Q2 = [A, C, D;
      C, B, E;
      D, E, F]

figure(3)
plotConic(Q2,'go');

P = intersectConics(Q1,Q2)
figure(1)
plot(P(1,:) ./ P(3,:) , P(2,:) ./ P(3,:), 'ro');
hold off;
