function sem = sem_weight(x, w)
    n = length(w);
    xWbar = sum(x .* w)/sum(w);
    wbar = mean(w);
    sem = n/((n-1)*sum(w)^2)*(sum((w.*x - wbar*xWbar).^2) - ...
            2*xWbar.*sum((w - wbar) .* (w.*x - wbar.*xWbar)) + xWbar.^2.*sum((w-wbar).^2));
    sem = sqrt(sem);
end