function set_range(rng)

set(0, 'currentfigure', 1);
h = get(1,'Children');

if isunix
    pre = findobj('UserData','pre_ax');
    set(pre(1),'clim',rng); set(pre(2),'clim',rng);
    pre = findobj('UserData','copy_pre_ax');
    set(pre(1),'clim',rng); set(pre(2),'clim',rng);
    push = findobj('UserData','push_ax');
    set(push(1),'clim',rng); set(push(2),'clim',rng);
    push = findobj('UserData','copy_push_ax');
    set(push(1),'clim',rng); set(push(2),'clim',rng);
    cb = findobj('UserData','disp_cb');
    set(cb,'Limits',rng);
    res = findobj('UserData','res_ax');
    set(res,'ylim',rng);
elseif ispc
    pre = findobj('UserData','pre_ax');
    set(pre,'clim',rng); set(h(5),'clim',rng);
    push = findobj('UserData','push_ax');
    set(push,'clim',rng); set(h(3),'clim',rng);
    set(h(1),'ylim',rng);
end

end
