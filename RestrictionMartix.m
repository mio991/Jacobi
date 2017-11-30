function R=RestrictionMartix(n)
  R=[];
  for i = 1:1:2^(n-1)-1
    z = zeros(1,2^n-1);
    z(2*i-1) = 1/4;
    z(2*i) = 1/2;
    z(2*i+1) = 1/4;
    R=[R;z];
  endfor
endfunction