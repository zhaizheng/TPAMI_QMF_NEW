t = -3:0.05:3;
X = [t;sin(t)]+0.02*randn(2,length(t));
k = 20;
d = 1;
Output = Spherelet(X, X, k, d);
subplot(1,2,1)
plot(X(1,:),X(2,:),'*')
subplot(1,2,2)
plot(Output(1,:),Output(2,:),'*')
