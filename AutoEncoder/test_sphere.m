k = 8;
d = 1;
X1 = readNPY('Noisy.npy');
ind2  =  randperm(n_sample,n_interest);
target2 = X1(:,ind2);
Sphere2 = Spherelet(target2, X1, k, d);