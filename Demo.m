l=6;
h=1/(2^l);
x=(h:h:(2^l-1)*h)';
f=pi^2*sin(pi*x)
u0=zeros(2^l-1,1);
[u,A]=Helmholz(f,u0,3,3,eps)