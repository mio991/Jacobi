# Jacobi Methode
#
# Fuer das Numerik I Seminar.

function res = Jacobi(A, b, x0=zeros(size(b)), eps=0.00001, omega = 1, doTests=true)
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
  M = inv(B)*(B-omega*A)
  Nb = omega*inv(B)*b;
  
  # Calculate iteration depth if eps < 1 else depth = eps
  depth = 0;
  if eps < 1
    nor = 0;
    #Chose Norm
    if r
      nor = inf;
    elseif c
      nor = 1;      
    elseif s
      nor = "fro";      
    endif
    
    #TODO: x1 missing
    
    #Calculate needed iteration depth.
    q = norm(A, nor);
    depth = ceil (log(eps*(1-q)/norm(b, nor))/log(q));
  else
    depth = ceil (eps);
  endif
  
  # Iterate through Jacobi
  res = iterate(zeros(size(b), 1), M, Nb, depth);
endfunction

function x = iterate(x, M, Nb, stepsleft)
  x=M*x+Nb;
  stepsleft--;
  if(stepsleft > 0)
    x = iterate(x, M, Nb, stepsleft);
  endif
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