function q_out = qinver(q_in)
    q_out = qconj(q_in) ./ sum(qnorm(q_in))^2;
end