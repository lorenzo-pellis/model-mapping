function map_gen = map_gen_1to2(n,k)

% This function maps a single index in a group of size n to a pair (s,i)
% The inverse of this map is map_gen_2to1
% The counting order in a household of size 2 is:
% 1 <--> (0,0)
% 2 <--> (0,1)
% 3 <--> (0,2)
% 4 <--> (1,0)
% 5 <--> (1,1)
% No number is associated with (2,0), because it's a disease-free equilibrium
% 
% Last update: 12-05-2019 

max = (n+2)*(n+1)/2-1;

if k > max
    error(['Troubles: k too big! max = ',num2str(max)])
end

h = k;
for m=(n+1):-1:1
    if h>m
        h=h-m;
    else
        break;
    end
end
s=n+1-m;
i=h-1;
map_gen = [s i];


