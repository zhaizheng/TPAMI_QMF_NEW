figure('Position', [1, 1, 550, 400])
x = -3:0.1:3;
y = -3:0.1:3;
[X, Y] = meshgrid(x,y);
Z = exp(1-(X.^2+2*Y.^2));%*((X/5-Y/3).^2);
t = tiledlayout(2,3,"TileSpacing","compact");

nexttile
q = 1;
Z = exp(1-(X.^2+2*Y.^2)).^q;
mesh(X,Y,Z)
hold on
Y1 = zeros(size(x));
Z1 = exp(1-(x.^2+2*Y1.^2)).^q;
plot3(x,Y1,Z1,'r-','Linewidth',2)
hold on
Y2 = [-3:0.1:-1/sqrt(8*q),1/sqrt(8*q):0.1:3];
X2 = zeros(size(Y2));
Z2 = exp(1-(X2.^2+2*Y2.^2)).^q;
plot3(X2,Y2,Z2,'r-','Linewidth',2)
axis([-3 3 -3  3 0 2.8])
title('q=1')


nexttile
q = 1/3;
Z = exp(1-(X.^2+2*Y.^2)).^q;
mesh(X,Y,Z)
hold on
Y1 = zeros(size(x));
Z1 = exp(1-(x.^2+2*Y1.^2)).^q;
plot3(x,Y1,Z1,'r-','Linewidth',2)
hold on
Y2 = [-3:0.1:-1/sqrt(8*q),1/sqrt(8*q):0.1:3];
X2 = zeros(size(Y2));
Z2 = exp(1-(X2.^2+2*Y2.^2)).^q;
plot3(X2,Y2,Z2,'r-','Linewidth',2)
axis([-3 3 -3  3 0 2.8])
title('q=1/3')


nexttile
q = 1/30;
Z = exp(1-(X.^2+2*Y.^2)).^q;
mesh(X,Y,Z)
hold on
Y1 = zeros(size(x));
Z1 = exp(1-(x.^2+2*Y1.^2)).^q;
plot3(x,Y1,Z1,'r-','Linewidth',2)
hold on
Y2 = [-3:0.1:-1/sqrt(8*q),1/sqrt(8*q):0.1:3];
X2 = zeros(size(Y2));
Z2 = exp(1-(X2.^2+2*Y2.^2)).^q;
plot3(X2,Y2,Z2,'r-','Linewidth',2)
axis([-3 3 -3  3 0 2.8])
title('q=1/30')




nexttile
q = -1/16;
Z = -exp(1-(X.^2+2*Y.^2)).^q;
mesh(X,Y,Z)
hold on
x = -sqrt(-1/2/q):0.1:sqrt(-1/2/q);
Y1 = zeros(size(x));
Z1 = -(exp(1-(x.^2))).^(q);
plot3(x,Y1,Z1,'r-','Linewidth',2)
title('q=-1/16')

nexttile
q = -1/8;
Z = -exp(1-(X.^2+2*Y.^2)).^q;
mesh(X,Y,Z)
hold on
x = -sqrt(-1/2/q):0.1:sqrt(-1/2/q);
Y1 = zeros(size(x));
Z1 = -(exp(1-(x.^2))).^(q);
plot3(x,Y1,Z1,'r-','Linewidth',2)
title('q=-1/8')

nexttile
q = -1;
Z = -exp(1-(X.^2+2*Y.^2)).^q;
mesh(X,Y,Z)
hold on
x = -sqrt(-1/2/q):0.1:sqrt(-1/2/q);
Y1 = zeros(size(x));
Z1 = -(exp(1-(x.^2))).^(q);
plot3(x,Y1,Z1,'r-','Linewidth',2)
title('q=-1')

exportgraphics(t,'demo_simple.eps','Resolution',400)




% %%
% syms x y q;
% z = exp(1-x^2-2*y^2);
% g = gradient(z);
% H = jacobian(g);
% % TH = H/z-(1-q) *g*g'/;
% pretty(g/z)
% pretty(H)