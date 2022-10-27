function [x1,x2, x3, ThreeRoots] = SolveCubic(A, B, C)
  p = (3*B-(A*A))/3;
  q = (2*(A*A*A)-(9*A*B)+(27*C))/27;
  m = 2*sqrt(-p/3);
  theta = acos(3*q/(p*m))/3;
  
  R = (q*q/4)+(p*p*p/27);
  
  if(R<0)
    x1 = m*cos(theta);
    x2 = m*cos(theta + 4*pi/3);
    x3 = m*cos(theta + 2*pi/3);
    ThreeRoots = 1;
  else
    x1 = 0;
    x2 = 0;
    x3 = 0;
    ThreeRoots = 0;
  endif
  
  
endfunction
