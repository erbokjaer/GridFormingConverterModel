function y = fault_profile(t, p)

t_fault   = p.t_fault;
t_apply   = p.t_apply;
t_clear   = p.t_clear;
t_recover = p.t_recover;

y_pre  = p.y_prefault;
y_fault= p.y_fault;
y_post = p.y_postfault;

if t < t_fault

    y = y_pre;

elseif t <= t_fault + t_apply

    y = y_pre + (y_fault - y_pre) * (t - t_fault)/t_apply;

elseif t < t_clear

    y = y_fault;

elseif t <= t_clear + t_recover

    y = y_fault + (y_post - y_fault) * (t - t_clear)/t_recover;

else

    y = y_post;

end

end