function uFein = twogrid(fFein, AFein, uFein, AGrob, pre = 5, post = 5, tol=0.0001, maxit=10000, omega=1/2)
  # Prepwork
  fGrob = LinearRestriction(fFein);
  figure(1);
  plot(uFein);
  # Calculations
  uFein = Jacobi(AFein, fFein, uFein, eps, pre, omega, false);
  figure(2);
  plot(uFein);
  uGrob = LinearRestriction(uFein);
  uGrob = Jacobi(AGrob, fGrob, uGrob, tol, maxit, omega, false);
  uFein = LinearProlongation(uGrob);
  figure(3);
  plot(uFein);
  uFein = Jacobi(AFein, fFein, uFein, eps, post, omega, false);
  figure(4);
  plot(uFein);
endfunction