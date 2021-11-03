function sample = tnrnd(m, s, range, n)
    untruncated = makedist('Normal', m, s);
    truncated = truncate(untruncated, range(1), range(2));
    sample = random(truncated, 1, n);
end
