function [x1,x2, x3] = SolveCubic(A, B, C)
  B = - B;
  Q = ((3*B)-(A*A))/9;
  R = ((9*A*B)-(27*C)-(2*A*A*A))/54;
  D = ((Q*Q*Q)+(R*R));
  
  theta = acos(R/(sqrt(-(Q*Q*Q))));
  
  
  if(D<0)
    x1 = 2*sqrt(-Q)*cos(theta/3)-(A/3)
    x2 = 2*sqrt(-Q)*cos(theta/3+ (2*pi/3))-(A/3)
    x3 = 2*sqrt(-Q)*cos(theta/3+ (4*pi/3))-(A/3)
    return;
  else
    x1 = 0;
  endif
endfunction
