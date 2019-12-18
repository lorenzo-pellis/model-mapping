function esc = esc_prob_distr(w,PI_single,Rh,Rg)

max_size = length(PI_single);
esc_p = 0;
for n = 1:max_size
    if PI_single(n) ~= 0
        esc_p = esc_p + PI_single(n) * pgfBall(exp(-Rg*w),n,1,Rh,0,1,1); % Careful: here eta = 1, but I don't want to generalise it...
    end
end
esc = esc_p;

