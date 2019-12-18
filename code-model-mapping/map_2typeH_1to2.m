function pair = map_2typeH_1to2(h,w)

% This is a hand-made function to map a single index into a pair of indices
% to loop through all possible household compositions in a single loop. I
% only later discovered Matlab has a built-in function that does the same
% job (called ind2sub), but it seems this approach is faster (at least for
% the small numbers I work with), so I'm just leaving it as it is.
% 
% Input:
%   - h = [h_a,h_c]: household composition, which can contain 0s, so I'll
%   have to work by increasing all indices by 1
%   - w: single index
% 
% Output: pair (a,c) corresponding to the single index w, with a from 0 to 
% h_a and c from 0 to h_c. In counting, I first loop through c for each a:
%
% (0,0) <--> 1
% (0,1) <--> 2
%  ...
% (0,h_c) <--> h_c+1
% (1,0) <--> 1*(h_c+1) + 1
% (1,1) <--> 1*(h_c+1) + 2
%  ...
% (1,h_c) <--> 1*(h_c+1) + (h_c+1) = 2*(h_c+1)
%  ...
% (h_a,h_c) <--> h_a*(h_c+1) + (h_c+1) = (h_a+1)*(h_c+1)
%
% Author: Lorenzo Pellis
% Last update: 20/05/2019 

max = (h(1)+1)*(h(2)+1); % Maximum value of the single index
if w > max
    error('Error: too large index!');
end
% I could use sume form of integer division and rest, but I'll just loop
a = 0;
l = w;
jump = h(2) + 1;
l = l - jump;
while l>0
    l = l - jump;
    a = a + 1;
end
c = l + jump - 1;
pair = [ a, c ];
