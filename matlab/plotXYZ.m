infile = '../mcraman/skinvessel_XYZ.txt';
data = load(infile, '-ascii');
x = data(:,1); % column 1
y = data(:,2); % column 2
z = data(:,3); % column 3
myColor = data(:,4); % column 4

figure
scatter3(x,y,z);
