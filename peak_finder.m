Prob = [1;2;3;4;5;4;3;2;1;2;3;4;5;6;7;8;9;8;7;6;5;4;3;2;1];
[arlen, arwid] = size(Prob);    %Get array length of 1D array
maxima = [];                    %Instantiate maxima list

for i=2:arlen-2                 %Iterate over indices with indices = index-1 and index+1 present
   if Prob(i) > Prob(i-1) && Prob(i) > Prob(i+1)
       maxima = [maxima; i];    %If value at index i is greater than adjacent values, add index to maxima list
   end
end
maxima_vals = Prob(maxima);     %Values at maxima

max1 = maxima(maxima_vals==max(maxima_vals));     %Index of largest maxima
max2 = maxima(maxima_vals==max(maxima_vals(maxima_vals<max(maxima_vals))));   %Index of 2nd largest maxima