
% Convert each diameter value from binary to decimal

function R = base2todec(P,nb)

[n,N] = size(P);

s = 1;
R = nan(n/nb,N);

% Convert each diameter value from binary to decimal
while(s <= n)
    for j = 1:n/nb
        R(j,:) = 2.^(nb-1:-1:0)*P(s:s+nb-1,:) + 1;
        s = s + nb;   
    end
end


