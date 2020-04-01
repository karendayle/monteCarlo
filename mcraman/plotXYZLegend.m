% Colors:
blue =    [0.0000, 0.4470, 0.7410];
rust =    [0.8500, 0.3250, 0.0980];
gold =    [0.9290, 0.6940, 0.1250];
purple =  [0.4940, 0.1840, 0.5560];
green =   [0.4660, 0.6740, 0.1880];
ciel =    [0.3010, 0.7450, 0.9330];
cherry =  [0.6350, 0.0780, 0.1840];
red =     [1.0, 0.0, 0.0];
black =   [0.0, 0.0, 0.0];

myTextFont = 30;
    
figure

x = 0.2;
y = 0.95;
z = 1.0;
deltaY = 0.15;
text(x, y, z, 'no change', 'Color', black, 'FontSize', myTextFont);
text(x, y, z, '_________', 'Color', black, 'FontSize', myTextFont);
y = y - deltaY;
text(x, y, z, '1082 cm-1', 'Color', red, 'FontSize', myTextFont);
text(x, y, z, '_________', 'Color', red, 'FontSize', myTextFont);
y = y - deltaY;
text(x, y, z, '1430 cm-1', 'Color', green, 'FontSize', myTextFont);
text(x, y, z, '_________', 'Color', green, 'FontSize', myTextFont);
y = y - deltaY;
text(x, y, z, '1584 cm-1', 'Color', blue, 'FontSize', myTextFont);
text(x, y, z, '_________', 'Color', blue, 'FontSize', myTextFont);
y = y - deltaY;
text(x, y, z, '1702 cm-1', 'Color', purple, 'FontSize', myTextFont);
text(x, y, z, '_________', 'Color', purple, 'FontSize', myTextFont);
y = y - deltaY;
text(x, y, z, 'other values', 'Color', rust, 'FontSize', myTextFont);
text(x, y, z, '__________', 'Color', rust, 'FontSize', myTextFont);
y = y - deltaY;

    