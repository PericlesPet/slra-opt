function P = stiefSVD_proj(R)

    [U, ~, V] = svd(R);
    P = U*eye(size(R))*V';

end