function y = ramp_signal(t, t_start, t_dur, y0, y1)

if t_dur == 0

    if t < t_start
        y = y0;
    else
        y = y1;
    end
    return

end

if t < t_start
    y = y0;

elseif t <= (t_start + t_dur)
    y = y0 + (y1 - y0) * (t - t_start) / t_dur;

else
    y = y1;

end

end