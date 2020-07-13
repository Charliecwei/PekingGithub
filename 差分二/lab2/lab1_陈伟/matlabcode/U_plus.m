function u = U_plus(u)
%% The positive part of u
    a = zeros(size(u));
    u = max([u';a'])';
end
