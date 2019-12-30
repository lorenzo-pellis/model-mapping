function av_gen_chain = create_chain_1type_check(PI,Rh,eta,str_den)

% This function computes the average number of cases in each generation
% (called chain) for a households model with households of variable size.
% The computation is done in 2 (or 3) ways: first, using a matrix
% representation of the probabilities of each possible Reed-Frost epidemic
% ("matrix_avRFchain_den"), but I don't comment on it, as it is not the 
% most efficient way; the other method is a recursive relation for the mean
% numbers, which is faster and is what I do in the function
% "create_chain_1type" I use in the main code. (The third version I try is
% the same as this second one, but passing the matrix of binomial
% coefficients, rather than declaring it a global variable: it is a tiny
% bit faster, but I'd have to change all similar codes - this is what I use
% now below).

max_gen = length(PI);
av_gen_chain_temp = zeros(1,max_gen);
% tic;
for n = 1:max_gen
    if PI(n) ~= 0
        av_gen_chain_temp(1:n) = av_gen_chain_temp(1:n) + PI(n) * matrix_avRFchain_den(n,1,Rh,0,1,1,'n-1'); % Careful: here eta = 1, and don't want to generalise
    end
end
av_gen_chain = av_gen_chain_temp;
% toc;


% av_gen_chain_temp = zeros(1,max_gen);
% % tic;
% for n = 1:max_gen
%     if PI(n) ~= 0
%         if strcmp(str_den,'n-1')
%             lambda = Rh/(n-1)^eta; % 1-to-1 escaping probability
%         else
%             lambda = Rh/n^eta;
%         end
%         av_gen_chain_temp(1:n) = av_gen_chain_temp(1:n) + PI(n) * create_avRFchain_3dim(n,1,lambda);
%     end
% end
% av_gen_chain2 = av_gen_chain_temp;
% % toc;
% % disp( [ av_gen_chain, av_gen_chain2 ] )


global bincoeffmat;

av_gen_chain_temp = zeros(1,max_gen);
% tic;
for n = 1:max_gen
    if PI(n) ~= 0
        if strcmp(str_den,'n-1')
            lambda = Rh/(n-1)^eta; % 1-to-1 escaping probability
        else
            lambda = Rh/n^eta;
        end
        av_gen_chain_temp(1:n)=av_gen_chain_temp(1:n) + PI(n) * create_avRFchain_3dim_pass(n,1,lambda,bincoeffmat);
    end
end
av_gen_chain3 = av_gen_chain_temp;
% toc;
% disp([ av_gen_chain; av_gen_chain3 ])

