% ref https://www.mathworks.com/matlabcentral/answers/486172-how-can-i-plot-a-filled-cylinder-with-a-specific-height
% Values for 2 cylinders
r = [0.05, 1];             % radii of each cyl
cnt = [1,0.2; 3,5];       % [x,y] center of each cyl
height = [1,8];         % height of each cyl
color = [1 0 0; 0 .5 0];% color of each cyl
nSides = 100;           % number of "sides" of the cyl
% Create figure
figure()
hold on
% Loop through each cylinder
for i = 1:1 % numel(r)
    plotCylinderWithCaps(r(i),cnt(i,:),height(i),nSides,color(i,:));
end
view(3)
grid on
function [h1, h2, h3] = plotCylinderWithCaps(r,cnt,height,nSides,color)
[X,Y,Z] = cylinder(r,nSides);
capColor = [0., 0., 0.]; % black
X = X + cnt(1); 
Y = Y + cnt(2); 
%Z = Z * height; 
Z = Z* height-1*height/2; 
h1 = surf(X,Z,Y,'facecolor',color,'LineStyle','none');
h2 = fill3(X(1,:),Z(1,:),Y(1,:),capColor);
h3 = fill3(X(2,:),Z(2,:),Y(2,:),capColor);
xlabel('x [cm]');
ylabel('y [cm]');
zlabel('z [cm]');
end  %only needed if this is within a script