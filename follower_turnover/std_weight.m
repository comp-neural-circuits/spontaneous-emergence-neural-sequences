function std = std_weight(x, w)
        xWbar = sum(x .* w)/sum(w);
        std = sqrt(sum(w.*(x - xWbar).^2)/sum(w));
end