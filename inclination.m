function theta = inclination(alpha,beta) % theta - inclination wrt + x axis
if beta <0 && beta >-pi/2
    theta = alpha - beta;
elseif beta <-pi/2 && beta >-pi
    if alpha + beta < pi
        theta = alpha - beta;
    else
        theta = alpha - beta - 2*pi;
    end
elseif beta >pi/2 && beta< pi
    if beta - alpha < pi
        theta = alpha - beta;
    else
        theta = alpha - beta + 2*pi;
    end
else
    theta = alpha - beta;
end
end

