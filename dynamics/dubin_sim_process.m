function x = dubin_sim_process(x, ts)
    global am
    if abs(x(3))<ts*am/2
        x(3)=0;
    end
end