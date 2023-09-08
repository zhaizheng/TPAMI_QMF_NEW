subplot(1,3,1)
x = sort(rand(1,10))-0.5;
z = sort(rand(1,100)-0.5);
y = x.^2+0.03*(rand(1,10)-0.5);
zy = z.^2;
zy_r = z.^2+0.03*(rand(1,100)-0.5);
plot(x,y,'-*');
subplot(1,3,2)
plot(z,zy_r,'o');
hold on
plot(z,zy,'-')


% A = rand(2,5);
% plot(A(1,:),A(2,:),'')