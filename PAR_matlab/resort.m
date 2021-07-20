function [Qn, Nfloat] = resort(Q, a, index)
n = size(Q,1);
Qn = zeros(n, n);
Nfloat = zeros(n, 1);
for i=1:n
    pos1 = find(index==i);
    Qn(pos1, pos1) = Q(i,i);
    Nfloat(pos1) = a(i);
    for j=1:n
        if j==i
            continue
        end
        pos2 = find(index==j);
        Qn(pos1,pos2) = Q(i, j);
        Qn(pos2,pos1) = Q(j, i);
    end
end