function u = U_plus_delta(u)
%% The positive part of u for parameter delta
    a = zeros(size(u));
    u = max_delta(u,a);
end
