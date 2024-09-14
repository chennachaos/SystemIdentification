%
% Standard regression technique for identifying
% Lorenz system from time series data
%
% Author: Dr Chennakesava Kadapa
%
% Date: 12-Sep-2024
%


clear all;    clc;
dt  =  0.01;
y0 = [-8, 8, 27];

myfun = @(t,y) [-10*y(1)+10*y(2);
                28.0*y(1)-y(1)*y(3)-y(2);
                y(1)*y(2)-8.0*y(3)/3.0];

[t,data] = ode45(myfun, [0.0:dt:40.0], y0);

x = data(:,1);
y = data(:,2);
z = data(:,3);

% find the derivative terms
xd = -10.0*x + 10.0*y;
yd = 28.0*x - x.*z - y;
zd = x.*y - 8.0*z./3.0;

% Construct matrix X
N = size(t,1);
A = zeros(3*N,30);
vvec = zeros(3*N,1);

for i=1:N
%for i=3:N
  ind1 = 3*(i-1)+1;
  ind2 = ind1+1;
  ind3 = ind1+2;

  A(ind1, 1)  = 1.0;
  A(ind1, 2)  = x(i);
  A(ind1, 3)  = y(i);
  A(ind1, 4)  = z(i);
  A(ind1, 5)  = x(i)*x(i);
  A(ind1, 6)  = y(i)*y(i);
  A(ind1, 7)  = z(i)*z(i);
  A(ind1, 8)  = x(i)*y(i);
  A(ind1, 9)  = y(i)*z(i);
  A(ind1, 10) = z(i)*x(i);

  A(ind2, 11:20) = A(ind1, 1:10);
  A(ind3, 21:30) = A(ind1, 1:10);

  vvec(ind1) = xd(i);
  vvec(ind2) = yd(i);
  vvec(ind3) = zd(i);
end

format long g

coeffs = inv(A'*A)*A'*vvec;
coeffs = reshape(coeffs,10,3)'


# testing with the coefficients computed
#
#y0 = data(end,:);

%[t1,data1] = ode45(myfun, [0:dt:40.0], y0);
t1 = t;
data1 = data;

##myfun2 = @(t,y) [coeffs(1,1)+coeffs(1,2)*y(1)+coeffs(1,3)*y(2)+coeffs(1,4)*y(3)+coeffs(1,5)*y(1)*y(1)+coeffs(1,6)*y(2)*y(2)+coeffs(1,7)*y(3)*y(3)+coeffs(1,8)*y(1)*y(2)+coeffs(1,9)*y(2)*y(3)+coeffs(1,10)*y(1)*y(3);
##                 coeffs(2,1)+coeffs(2,2)*y(1)+coeffs(2,3)*y(2)+coeffs(2,4)*y(3)+coeffs(2,5)*y(1)*y(1)+coeffs(2,6)*y(2)*y(2)+coeffs(2,7)*y(3)*y(3)+coeffs(2,8)*y(1)*y(2)+coeffs(2,9)*y(2)*y(3)+coeffs(2,10)*y(1)*y(3);
##                 coeffs(3,1)+coeffs(3,2)*y(1)+coeffs(3,3)*y(2)+coeffs(3,4)*y(3)+coeffs(3,5)*y(1)*y(1)+coeffs(3,6)*y(2)*y(2)+coeffs(3,7)*y(3)*y(3)+coeffs(3,8)*y(1)*y(2)+coeffs(3,9)*y(2)*y(3)+coeffs(3,10)*y(1)*y(3)];

[t2,data2] = ode45(@myfun2, [0:dt:40.0], y0);


%set(gcf, 'PaperUnits', 'inches');
%set(gcf, 'PaperSize', [10 6]);

subplot(3,1,1)
plot(t1,data1(:,1), 'b', "linewidth",2)
hold on
plot(t2,data2(:,1), 'k--', "linewidth",2)
ylabel("x")
xticks([])
legend('Truth','Model','Location','northwest')
title('With Octave')
pbaspect([6 1 1])

subplot(3,1,2)
plot(t1,data1(:,2), 'b', "linewidth",2)
hold on
plot(t2,data2(:,2), 'k--', "linewidth",2)
ylabel("y")
xticks([])
pbaspect([6 1 1])

subplot(3,1,3)
plot(t1,data1(:,3), 'b', "linewidth",2)
hold on
plot(t2,data2(:,3), 'k--', "linewidth",2)
xlabel("t")
ylabel("z")
xticks(0:5:40)
pbaspect([6 1 1])

saveas(gcf,'plot.png')

























