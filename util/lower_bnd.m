function z = lower_bnd(c, C0, tol, slope)
    % c = a matrix holding (x0, x1, x2, ..., xN)
    % c(1:end-1) --> x0 to x(N-1) are parameter inputs
    % c(end) --> xN is the output
    if strcmp(c, 'init')
    else
        [ ~ , cols] = size(C0);           
        z = max( (C0(end,:) - tol) - slope*sqrt(sum((c*ones(1,cols) - C0(1:end-1,:)).^2)) );
    end
end