function map_gen = map_gen_2to1(n,s,i)

% This function maps a pair (s,i) in a group of size n into a single index
% The inverse of this map is map_gen_1to2
% The counting order in a household of size 2 is:
% 1 <--> (0,0)
% 2 <--> (0,1)
% 3 <--> (0,2)
% 4 <--> (1,0)
% 5 <--> (1,1)
% No number is associated with (2,0), because it's a disease-free equilibrium
% 
% Last update: 12-05-2019 

h = 0;
% max = (n+2)*(n+1)/2-1

if (s+i) > n
    error(['Troubles: invalid pair of susc-inf! n = ',num2str(n)]);
end
if s == n
    error(['i = ',num2str(i),'...  Troubles: everybody is susceptible!']);
end

for m=(n+1):-1:(n+2-s)
    h = h+m;
end
k=h+i+1;
map_gen = k;


