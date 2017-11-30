function x = LinearRestriction(x)
  [a,_] = size(x);
  n = log2(a+1);
  
  R=RestrictionMartix(n);  
  x=R*x;
endfunction