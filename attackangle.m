function alpha = attackangle(beta,theta)
if beta <0 && beta >-pi/2
    %theta = alpha - beta;
    alpha = beta + theta;
elseif beta <-pi/2 && beta >-pi
    if theta > 0
        %theta = alpha - beta;
        alpha = beta + theta;
    else
        %theta = alpha - beta - 2*pi;
        alpha = beta + theta + 2*pi;
    end
elseif beta >pi/2 && beta< pi
    if theta < 0
        %theta = alpha - beta;
        alpha = theta + beta;
    else
        %theta = alpha - beta + 2*pi;
        alpha = beta + theta - 2*pi;
    end
else
    %theta = alpha - beta;
    alpha = beta + theta;
end
end