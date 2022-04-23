%% trajectory generation
clear all;
close all;
% generate robot1 local trajectory
pose_robot1 = [0 5 7 16 20;
               0 10 19 26 32;];
angle_robot1 = [0 pi/3 pi/2 pi/4 pi/6];
plot(pose_robot1(1,:),pose_robot1(2,:),'r')
hold on;
% generate robot1 global trajectory
global_pose_robot1 = pose_robot1;

% generate robot2 local trajectory
pose_robot2 = [0 3 10 17 22;
               0 11 21 26 36;];
angle_robot2 = [0 pi/4 pi/5 pi/6 pi/3];
% plot(pose_robot2(1,:),pose_robot2(2,:),'y')

% translation from robot2's frame to robot1's frame
translation = [1; 6];
true_angle = atan(translation(2)/translation(1));

% rotation from robot2's frame to robot1's frame
% Coordinate system 2 rotates pi / 8 clockwise relative to coordinate system 1;
rotation_angle = pi/15;
rotation_matrix = [cos(rotation_angle) -sin(rotation_angle);
                   sin(rotation_angle) cos(rotation_angle)];

% generate global_pose_robot2
global_pose_robot2 = zeros(2,5);
for n = 1:5
    global_pose_robot2(:,n) = rotation_matrix*pose_robot2(:,n);
end
% plot(global_pose_robot2(1,:),global_pose_robot2(2,:),'r')

for n = 1:5
    global_pose_robot2(:,n) = global_pose_robot2(:,n) + translation;
end
plot(global_pose_robot2(1,:),global_pose_robot2(2,:),'g')

% compute distance
distance = zeros(1,5);

for n = 1:5
    tmp_distance = global_pose_robot1(:,n) - global_pose_robot2(:,n);
    distance(1,n) = norm(tmp_distance);
end

%% 
A = zeros(4,7);
for n =  2:5
    A(n-1,1) = pose_robot1(2,n)*pose_robot2(2,n)+pose_robot1(1,n)*pose_robot2(1,n);
    A(n-1,2) = pose_robot1(2,n)*pose_robot2(1,n)-pose_robot1(1,n)*pose_robot2(2,n);
    A(n-1,3) = distance(1,1)*pose_robot1(1,n);
    A(n-1,4) = distance(1,1)*pose_robot1(2,n);
    A(n-1,5) = -distance(1,1)*pose_robot2(1,n);
    A(n-1,6) = -distance(1,1)*pose_robot2(2,n);
    A(n-1,7) = 0.5*(distance(1,n)*distance(1,n)-distance(1,1)*distance(1,1)-transpose(pose_robot2(:,n))*pose_robot2(:,n)-transpose(pose_robot1(:,n))*pose_robot1(:,n));
end

% compute nullspace of A
N = null(A);
r = N(:,1);
s = N(:,2);
t = N(:,3);

% compute e
e1 = r(1)-r(7)*t(1)/t(7);
e2 = s(1)-s(7)*t(1)/t(7);
e3 = t(1)/t(7);
e4 = r(2)-r(7)*t(2)/t(7);
e5 = s(2)-s(7)*t(2)/t(7);
e6 = t(2)/t(7);

e7 = r(3)-r(7)*t(3)/t(7);
e8 = s(3)-s(7)*t(3)/t(7);
e9 = t(3)/t(7);
e10 = r(4)-r(7)*t(4)/t(7);
e11 = s(4)-s(7)*t(4)/t(7);
e12 = t(4)/t(7);

e13 = r(5)-r(7)*t(5)/t(7);
e14 = s(5)-s(7)*t(5)/t(7);
e15 = t(5)/t(7);
e16 = r(6)-r(7)*t(6)/t(7);
e17 = s(6)-s(7)*t(6)/t(7);
e18 = t(6)/t(7);

% compute beta
Beta = zeros(5,5);
Beta(1,1) = e1*e1 + e4*e4;
Beta(1,2) = e2*e2 + e5*e5;
Beta(1,3) = 2*(e1*e2 + e4*e5);
Beta(1,4) = 2*(e1*e3 + e4*e6);
Beta(1,5) = 2*(e2*e3 + e5*e6);

Beta(2,1) = e7*e7 + e10*e10;
Beta(2,2) = e8*e8 + e11*e11;
Beta(2,3) = 2*(e7*e8 + e10*e11);
Beta(2,4) = 2*(e7*e9 + e10*e12);
Beta(2,5) = 2*(e8*e9 + e11*e12);

Beta(3,1) = e13*e13 + e16*e16;
Beta(3,2) = e14*e14 + e17*e17;
Beta(3,3) = 2*(e13*e14 + e16*e17);
Beta(3,4) = 2*(e13*e15 + e16*e18);
Beta(3,5) = 2*(e14*e15 + e17*e18);

Beta(4,1) = e1*e7 + e4*e10;
Beta(4,2) = e2*e8 + e5*e11;
Beta(4,3) = e1*e8 + e2*e7 + e4*e11 + e5*e10;
Beta(4,4) = e1*e9 + e3*e7 + e4*e12 + e6*e10 - e13;
Beta(4,5) = e2*e9 + e3*e8 + e5*e12 + e6*e11 - e14;

Beta(5,1) = e1*e10 - e4*e7; 
Beta(5,2) = e2*e11 - e5*e8;
Beta(5,3) = e1*e11 + e2*e10 - e4*e8 - e5*e7;
Beta(5,4) = e1*e12 + e3*e10 - e4*e9 - e6*e7 - e16;
Beta(5,5) = e2*e12 + e3*e11 - e5*e9 - e6*e8 - e17;

% compute F
F = zeros(5,1);
F(1) = 1 - e3*e3 - e6*e6;
F(2) = 1 - e9*e9 - e12*e12;
F(3) = 1 - e15*e15 - e18*e18;
F(4) = e15 - e3*e9 - e6*e12;
F(5) = e18 - e3*e12 + e6*e9;

% compute y
y = Beta\F;

%compute lamda
lamda1 = y(4);
lamda2 = y(5);
lamda3 = (1 - lamda1*r(7) - lamda2*s(7))/t(7);

%compute x
x = lamda1*r+lamda2*s+lamda3*t;

% compute result
theta = asin(x(4));
phi = asin(x(2));
res_trans = [distance(1)*cos(theta) distance(1)*sin(theta)];