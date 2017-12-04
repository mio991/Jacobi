function u = Helmholz(f, u0)
  
  if(size(f) ~= size(u0))
    u=0;
    print("size(f) ~= size(u0)!\n");
    return;
  endif
  
  [a,_] = size(f);
  div = log2(a+1);
  
  AFein = Discretisierung(2^div-1);
  AGrob = Discretisierung(2^(div-1)-1);
  
  u = twogrid(f, AFein, u0, AGrob);
endfunction

function A = Discretisierung(n)
  D = (2*(n+1)^2-1)*eye(n);
  nD = -(n+1)^2*eye(n-1);
  OnD=[zeros(n-1,1), nD; zeros(1,n)];
  UnD=transpose(OnD);
  
  A=D+OnD+UnD;
endfunction