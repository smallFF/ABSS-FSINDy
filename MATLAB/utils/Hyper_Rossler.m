function dy = Hyper_Rossler(t,y,a,b,c,d)

dy = [
-y(2) - y(3);
y(1) + a*y(2) + y(4);
b + y(1)*y(3);
-d*y(3) + c*y(4)
];

