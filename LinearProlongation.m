function x = LinearProlongation(x)
  [a,_] = size(x);
  n = log2(a+1);
  
  P=transpose(2*RestrictionMartix(n+1));
  x=P*x;
endfunction