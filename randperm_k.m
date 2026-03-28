function seq = randperm_k(n, k)
% RANDPERM_K Generate a sequence of length k using repeated randperm(n)
%
% If k <= n, this is just randperm(n, k).
% If k > n, it concatenates multiple randperm(n) calls to ensure
% each block of n samples fully covers 1:n.

    if k <= n
        seq = randperm(n, k);
        return;
    end

    % Number of full permutations needed
    num_full = floor(k / n);
    remainder = mod(k, n);

    % Preallocate for speed
    seq = zeros(1, k);

    idx = 1;

    % Add full randperm blocks
    for i = 1:num_full
        rp = randperm(n);
        seq(idx:idx+n-1) = rp;
        idx = idx + n;
    end

    % Add remainder
    if remainder > 0
        rp = randperm(n);
        seq(idx:end) = rp(1:remainder);
    end
end