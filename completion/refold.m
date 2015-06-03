function out = refold(matrix, mode, Dims)
p = Dims(1); q = Dims(2); r = Dims(3);
switch mode
    case 1
        % Here 'dim' specifies 'r'
        out = zeros(p, q, r);
        for i = 1:r
            out(:, :, i) = matrix(:, q*(i-1)+1:q*i);
        end
    case 2
        % Here 'dim' specifies 'r'
        out = zeros(p, q, r);
        for i = 1: r
            out(:, :, i) = matrix(:, p*(i-1)+1:p*i)';
        end
    case 3
        % Here 'dim' specifies 'q'
        out = zeros(p, q, r);
        for i = 1:r
            out(:, :, i) = reshape(matrix(i, :), p, q);
        end
    otherwise
        error('Mode cannot be bigger than 3.')
end