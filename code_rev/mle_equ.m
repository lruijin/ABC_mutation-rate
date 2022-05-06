function mle = mle_equ(z,x,p)
    mle = (1-2*p)*(z-x) - (1-p)*z^(1-2*p)+p;
end