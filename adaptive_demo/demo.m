subplot(1,3,1)
x = 0:0.06:5;
y = x.^(0.5).*sin(2*x);
z = x.^(0.5).*sin(2*x)+0.15*randn(size(x));
plot(x,y)
hold on
plot(x,z,'*')

s = -2:0.3:2.5;

[X,Y] = meshgrid(x,s);
n = length(x)*length(s);
Z = zeros(length(s),length(x));
for i = 1:length(s)
    for j = 1:length(x)
        Z(i,j) = KDE([x;z],0.2,[X(i,j);Y(i,j)]);
    end
end
%Z = reshape(ZZ,[length(s),length(x)]);
subplot(1,3,2)
surf(X,Y,Z)
subplot(1,3,3)
surf(X,Y,Z)
grid off
box off


function z = KDE(x,h,a)
    z = 0;
    for i = 1:size(x,2)
        z = z+ exp(-norm(a-x(:,i))^2/h^2);
    end
end