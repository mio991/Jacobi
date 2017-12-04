# Jacobi Methode
#
# Fuer das Numerik I Seminar.

function res = Jacobi(A, b, x0=zeros(size(b)), epsilon=0.00001, maxit=10000, omega = 1, doTests=true)
  res = 0;
  # Pretest for the Jacobi Methode
  
  # Size equals
  
  [a1 , a2] = size(A);
  [v, _] = size(b);
  
  if(a1~=a2 || a2~=v)
    printf("The matrix must be square, and have the same size as the vector!");
    return;
  endif
  
  # Diag non zeros
  if(sum(diag(A)~=0)~=size(A))
    printf("All diagonal elements must be ~=0!\n");
    return;
  endif
  
  if doTests 
    # Row-Sum-Criterea (Inf-Norm)
    r = RowSumTest(A);
    # Column-Sum Criterea (1-Norm)
    c = RowSumTest(transpose(A));
    # Square-Sum-criterea (Frobenius-Norm)
    s = SquareSumTest(A);
    
    if!(r || c || s)
      if !r
      printf("Row-Sum-Criterea Faild!\n")
    elseif !c
      printf("Column-Sum-Criterea Faild!\n")
    elseif !s
      printf("Square-Sum-Criterea Faild!\n")
    endif
    return
    endif
  endif
  
  # Prepare Jacobi Methode
  B = diag(diag(A));
  M = inv(B)*(B-omega*A);
  Nb = omega*inv(B)*b;
  
  nor = inf;
  if doTests
    #Chose Norm
    if r
      nor = inf;
    elseif c
      nor = 1;      
    elseif s
      nor = "fro";      
    endif
  endif
  
  #Calculate needed iteration depth.
  q = norm(A, nor);
  x1 = x=M*x0+Nb;
  depth = ceil (log(epsilon*(1-q)/norm(x1-x0, nor))/log(q));
  
  depth = abs(min(abs(depth), abs(maxit)));
  
  x = x1;
  
  # Iterate through Jacobi
  while depth > 0
    x=M*x+Nb;
    depth--;
  endwhile
  
  res = x;
endfunction

function res = RowSumTest(A)  
  strong = sum(A,2)<2*diag(A);
  weak = sum(A,2)<=2*diag(A);
  
  if sum(weak)~=size(weak)
    res = false;
    return;
  endif
  
  if sum(strong)==size(strong)
    res = true;
    return;
  endif
    
  if(sum(strong)==0)
    res = false;
    return;
  endif
  
  # Test for connectivity through algebraic connectivity.
  libNetwork = "./octave-networks-toolbox/";
  addpath (libNetwork);
  
  res = algebraicConnectivity(A) > 0;
  
  rmpath(libNetwork);
endfunction

function res = SquareSumTest(A)
  A = A.^2;
  D = diag(A);
  A = A - diag(D);
  S = sum(A , 2) ./ D;
  
  res = sum(S) < 1  ;
endfunction