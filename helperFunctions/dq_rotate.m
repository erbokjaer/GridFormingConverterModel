function x = dq_rotate(xd, xq, delta)
%DQ_ROTATE Rotate dq components into a new reference frame
%
% Inputs:
%   xd, xq  : 1xN (or Nx1) signals
%   delta   : 1xN (or Nx1) angle [rad]
%
% Outputs:
%   xdc, xqc : rotated components

    c = cos(delta);
    s = sin(delta);

    xdc =  c .* xd - s .* xq;
    xqc =  s .* xd + c .* xq;

    x = [xdc; xqc];
    
end