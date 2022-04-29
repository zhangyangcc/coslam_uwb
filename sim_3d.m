%% trajectory generation
clear all;
close all;
% generate robot1 local trajectory
pose_robot1 = [0 5 7 16 20 26;
               0 10 19 26 32 40;
               0 12 15 23 33 40];
plot3(pose_robot1(1,:),pose_robot1(2,:),pose_robot1(3,:),'r');
grid on;
hold on;
% generate robot1 global trajectory
global_pose_robot1 = pose_robot1;

% generate robot2 local trajectory
pose_robot2 = [0 3 10 17 22 31;
               0 11 21 26 36 47;
               0 5  13 23 29 42];
           
% translation from robot2's frame to robot1's frame
translation = [1; 6; 4;];

% rotation from robot2's frame to robot1's frame
% Coordinate system 2 rotates pi / 8 clockwise relative to coordinate system 1;
rotation_angle = pi/15;
rotation_matrix = [cos(rotation_angle) -sin(rotation_angle) 0;
                   sin(rotation_angle) cos(rotation_angle)  0;
                   0                   0                    1;];
% generate global_pose_robot2
global_pose_robot2 = zeros(3,6);
for n = 1:6
    global_pose_robot2(:,n) = rotation_matrix*pose_robot2(:,n) + translation;
end
plot3(global_pose_robot2(1,:),global_pose_robot2(2,:),global_pose_robot2(3,:),'g');

% compute distance
distance = zeros(1,6);

for n = 1:6
    tmp_distance = global_pose_robot1(:,n) - global_pose_robot2(:,n);
    distance(1,n) = norm(tmp_distance);
end
% compelete distance generation

%% solution for the problem
A = zeros(5,8);
for n =  2:6
    A(n-1,1) = -pose_robot1(1,n);
    A(n-1,2) = -pose_robot1(2,n);
    A(n-1,3) = pose_robot2(3,n) - pose_robot1(3,n);
    A(n-1,4) = -pose_robot1(1,n)*pose_robot2(1,n) - pose_robot1(2,n)*pose_robot2(2,n);
    A(n-1,5) = -pose_robot1(2,n)*pose_robot2(1,n) + pose_robot1(1,n)*pose_robot2(2,n);
    A(n-1,6) = pose_robot2(1,n);
    A(n-1,7) = pose_robot2(2,n);
    A(n-1,8) = -0.5*(distance(1,n)^2-distance(1,1)^2-transpose(pose_robot2(:,n))*pose_robot2(:,n)-transpose(pose_robot1(:,n))*pose_robot1(:,n)) - pose_robot1(3,n)*pose_robot2(3,n);
end

% compute nullspace of A
N = null(A);
r = N(:,1);
s = N(:,2);
t = N(:,3);

% compute Beta
Beta = zeros(5,5);
Beta(1,1) = r(4)^2 + r(5)^2 + r(8)^2*(t(4)^2+t(5)^2)/t(8)^2 - 2*r(8)*(r(4)*t(4)+r(5)*t(5))/t(8);
Beta(1,2) = s(4)^2 + s(5)^2 + s(8)^2*(t(4)^2+t(5)^2)/t(8)^2 - 2*s(8)*(s(4)*t(4)+s(5)*t(5))/t(8);
Beta(1,3) = 2*r(4)*s(4) + 2*r(5)*s(5) + 2*r(8)*s(8)*(t(4)^2+t(5)^2)/t(8)^2 - 2*(r(8)*s(4)*t(4)+r(4)*s(8)*t(4)+r(8)*s(5)*t(5)+r(5)*s(8)*t(5))/t(8);
Beta(1,4) = -2*r(8)*(t(4)^2+t(5)^2)/t(8)^2 + 2*(r(4)*t(4)+r(5)*t(5))/t(8);
Beta(1,5) = -2*s(8)*(t(4)^2+t(5)^2)/t(8)^2 + 2*(s(4)*t(4)+s(5)*t(5))/t(8);

Beta(2,1) = r(1)^2 + r(2)^2 + r(3)^2 + r(8)^2*(t(1)^2+t(2)^2+t(3)^2)/t(8)^2 - 2*r(8)*(r(1)*t(1)+r(2)*t(2)+r(3)*t(3))/t(8);
Beta(2,2) = s(1)^2 + s(2)^2 + s(3)^2 + s(8)^2*(t(1)^2+t(2)^2+t(3)^2)/t(8)^2 - 2*s(8)*(s(1)*t(1)+s(2)*t(2)+s(3)*t(3))/t(8);
Beta(2,3) = 2*r(1)*s(1) + 2*r(2)*s(2) + 2*r(3)*s(3) + 2*r(8)*s(8)*(t(1)^2+t(2)^2+t(3)^2)/t(8)^2 - 2*(r(8)*s(1)*t(1)+r(1)*s(8)*t(1)+r(8)*s(2)*t(2)+r(2)*s(8)*t(2)+r(8)*s(3)*t(3)+r(3)*s(8)*t(3))/t(8);
Beta(2,4) = -2*r(8)*(t(1)^2+t(2)^2+t(3)^2)/t(8)^2 + 2*(r(1)*t(1) + r(2)*t(2) + r(3)*t(3))/t(8);
Beta(2,5) = -2*s(8)*(t(1)^2+t(2)^2+t(3)^2)/t(8)^2 + 2*(s(1)*t(1) + s(2)*t(2) + s(3)*t(3))/t(8);

Beta(3,1) = r(1)*r(4) + r(2)*r(5) + r(8)^2*(t(1)*t(4) + t(2)*t(5))/t(8)^2 - r(8)*(r(4)*t(1)+r(5)*t(2)+r(1)*t(4)+r(2)*t(5))/t(8);
Beta(3,2) = s(1)*s(4) + s(2)*s(5) + s(8)^2*(t(1)*t(4) + t(2)*t(5))/t(8)^2 - s(8)*(s(4)*t(1)+s(5)*t(2)+s(1)*t(4)+s(2)*t(5))/t(8);
Beta(3,3) = r(4)*s(1) + r(5)*s(2) + r(1)*s(4) + r(2)*s(5) + 2*r(8)*s(8)*(t(1)*t(4) + t(2)*t(5))/t(8)^2 - (r(8)*s(4)*t(1) + r(4)*s(8)*t(1) + r(8)*s(5)*t(2) + r(5)*s(8)*t(2) + r(8)*s(1)*t(4) + r(1)*s(8)*t(4) + r(8)*s(2)*t(5) + r(2)*s(8)*t(5))/t(8);
Beta(3,4) = -r(6) - 2*r(8)*(t(1)*t(4)+t(2)*t(5))/t(8)^2 + (r(4)*t(1) + r(5)*t(2) + r(1)*t(4) + r(2)*t(5) + r(8)*t(6))/t(8);
Beta(3,5) = -s(6) - 2*s(8)*(t(1)*t(4)+t(2)*t(5))/t(8)^2 + (s(4)*t(1) + s(5)*t(2) + s(1)*t(4) + s(2)*t(5) + s(8)*t(6))/t(8);

Beta(4,1) = r(2)*r(4) - r(1)*r(5) + r(8)^2*(t(2)*t(4) - t(1)*t(5))/t(8)^2 + r(8)*(r(5)*t(1)-r(4)*t(2)-r(2)*t(4)+r(1)*t(5))/t(8);
Beta(4,2) = s(2)*s(4) - s(1)*s(5) + s(8)^2*(t(2)*t(4) - t(1)*t(5))/t(8)^2 + s(8)*(s(5)*t(1)-s(4)*t(2)-s(2)*t(4)+s(1)*t(5))/t(8);
Beta(4,3) = -r(5)*s(1) + r(4)*s(2) + r(2)*s(4) - r(1)*s(5) + 2*r(8)*s(8)*(t(2)*t(4) - t(1)*t(5))/t(8)^2 + (r(8)*s(5)*t(1) + r(5)*s(8)*t(1) - r(8)*s(4)*t(2) - r(4)*s(8)*t(2) - r(8)*s(2)*t(4) - r(2)*s(8)*t(4) + r(8)*s(1)*t(5) + r(1)*s(8)*t(5))/t(8);
Beta(4,4) = -r(7) - 2*r(8)*(t(2)*t(4)-t(1)*t(5))/t(8)^2 + (-r(5)*t(1) + r(4)*t(2) + r(2)*t(4) - r(1)*t(5) + r(8)*t(7))/t(8);
Beta(4,5) = -s(7) - 2*s(8)*(t(2)*t(4)-t(1)*t(5))/t(8)^2 + (-s(5)*t(1) + s(4)*t(2) + s(2)*t(4) - s(1)*t(5) + s(8)*t(7))/t(8);

Beta(5,1) = r(3)^2 + r(6)^2 + r(7)^2 + r(8)^2*(t(3)^2+t(6)^2+t(7)^2)/t(8)^2 - 2*(r(3)*r(8)*t(3)+r(6)*r(8)*t(6)+r(7)*r(8)*t(7))/t(8);
Beta(5,2) = s(3)^2 + s(6)^2 + s(7)^2 + s(8)^2*(t(3)^2+t(6)^2+t(7)^2)/t(8)^2 - 2*(s(3)*s(8)*t(3)+s(6)*s(8)*t(6)+s(7)*s(8)*t(7))/t(8);
Beta(5,3) = 2*r(3)*s(3) + 2*r(6)*s(6) + 2*r(7)*s(7) + 2*r(8)*s(8)*(t(3)^2+t(6)^2+t(7)^2)/t(8)^2 - 2*(r(8)*s(3)*t(3)+r(3)*s(8)*t(3)+r(8)*s(6)*t(6)+r(6)*s(8)*t(6)+r(8)*s(7)*t(7)+r(7)*s(8)*t(7))/t(8);
Beta(5,4) = -2*r(8)*(t(3)^2+t(6)^2+t(7)^2)/t(8)^2 + 2*(r(3)*t(3) + r(6)*t(6) + r(7)*t(7))/t(8);
Beta(5,5) = -2*s(8)*(t(3)^2+t(6)^2+t(7)^2)/t(8)^2 + 2*(s(3)*t(3) + s(6)*t(6) + s(7)*t(7))/t(8);

% compute F
F = zeros(5,1);
F(1) = 1 - (t(4)^2+t(5)^2)/t(8)^2;
F(2) = distance(1)^2 - (t(1)^2+t(2)^2+t(3)^2)/t(8)^2;
F(3) = -(t(1)*t(4)+t(2)*t(5))/t(8)^2+t(6)/t(8);
F(4) = -(t(2)*t(4)-t(1)*t(5))/t(8)^2+t(7)/t(8);
F(5) = distance(1)^2 - (t(3)^2+t(6)^2+t(7)^2)/t(8)^2;

% compute y
y = Beta\F;

%compute lamda
lamda1 = y(4);
lamda2 = y(5);
lamda3 = (1 - lamda1*r(8) - lamda2*s(8))/t(8);

%compute x
x = lamda1*r+lamda2*s+lamda3*t;

% compute result

phi = asin(x(5));
res_trans = x(1:3);