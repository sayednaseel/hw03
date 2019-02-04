function xerror = rimlesswheelerror(x, parms)
    if nargin < 2 % use default parameters
        parms = struct('alpha', 0.3, 'rgyr', 0, 'gamma', 0.04, 'tmax', 5);
    end
    
    [xnext, te, xs, ts, energies] = rimlesswheelstep(x, parms);
    xerror = xnext - x;
end