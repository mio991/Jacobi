function [u, AFein]= Helmholz(f, u0,pre=5,post=5,tol=0.0001,maxit=10000,omega=1/2)
  
  if(size(f) ~= size(u0))
    u=0;
    print("size(f) ~= size(u0)!\n");
    return;
  endif
  
  [a,_] = size(f);
  div = log2(a+1);
  
  AFein = Discretisierung(2^div-1);
  AGrob = Discretisierung(2^(div-1)-1);
  
  u = twogrid(f, AFein, u0, AGrob,pre,post,tol,maxit,omega);
endfunction

function A = Discretisierung(n)
  D = (2*(n+1)^2-1)*eye(n);
  nD = -(n+1)^2*eye(n-1);
  OnD=[zeros(n-1,1), nD; zeros(1,n)];
  UnD=transpose(OnD);
  
  A=D+OnD+UnD;
endfunction