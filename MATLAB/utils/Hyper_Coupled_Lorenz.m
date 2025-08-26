function dy = Hyper_Coupled_Lorenz(t,y,K,R)

dy = [
-y(1) + y(2) + K*y(4);
-y(1)*y(3);
-R + y(1)*y(2);
K*y(1) - y(4) + y(5);
-y(4)*y(6);
-R + y(4)*y(5);
];

