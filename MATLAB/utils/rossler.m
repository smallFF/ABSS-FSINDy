function dy = rossler(t,y,alpha,beta)
dy = [
-(y(2)+y(3));
y(1)+alpha*y(2);
alpha+y(3)*(y(1)-beta);
];

