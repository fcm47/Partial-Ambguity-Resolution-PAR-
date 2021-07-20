function [Qb1, Qab1, xb1] = deleteQ_Qab_xb(Qb, Qab, xb, list)
for i=length(list):-1:1
    ind = list(i);
    xb = [xb(1:ind-1); xb(ind+1:end)];
    Qab = [Qab(:, 1:ind-1), Qab(:, ind+1:end)];
    Qb = [Qb(:, 1:ind-1), Qb(:, ind+1:end)];
    Qb = [Qb(1:ind-1, :); Qb(ind+1:end, :)];
end

Qb1 = Qb;
Qab1 = Qab;
xb1 = xb;