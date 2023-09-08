 x = linspace(-1,1,100);
 y = linspace(-1,1,100);
 [X,Y] = meshgrid(x,y);
 Z = X.^2+3*Y.^2-8.*X.*Y;
 surf(X,Y,Z)

%%
[U,L,V] = svd([0,1;1,0]);