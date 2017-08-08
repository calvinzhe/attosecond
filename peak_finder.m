Prob = [1;2;3;4;5;4;3;2;1;2;3;4;5;6;7;8;9;8;7;6;5;4;3;2;1];
i_max = find(Prob==max(Prob));
i_nextmax = 0;
temp_Prob = Prob;
i_loop = 0;
while (i_nextmax == 0) && (i_loop < 11)
    i_currmax_t = find(max(temp_Prob));
    i_nextmax_d = find(Prob==max(Prob(Prob<max(temp_Prob))));
    [len,wid]=size(i_nextmax_d);
    if len > 1
       i_nextmax_d = i_nextmax_d(1); 
    end
    if (Prob(i_nextmax_d) > Prob(i_nextmax_d+1)) && (Prob(i_nextmax_d) > Prob(i_nextmax_d-1))
        i_nextmax = i_nextmax_d;
    else
        [Prob(i_nextmax_d), Prob(i_nextmax_d+1), Prob(i_nextmax_d-1)]
        temp_Prob(i_nextmax_d) = 0;
        temp_Prob(i_currmax_t) = 0;
        temp_Prob
    end
    i_loop = i_loop + 1;
    figure, plot(temp_Prob)
end
i_max
i_nextmax
Prob(i_max)
Prob(i_nextmax)
i_loop
