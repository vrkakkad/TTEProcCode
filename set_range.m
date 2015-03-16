function set_range(rng)

set(0, 'currentfigure', 1);
h = get(1,'Children');

pre = findobj('UserData','pre_ax');
set(pre,'clim',rng); set(h(5),'clim',rng);
push = findobj('UserData','push_ax');
set(push,'clim',rng); set(h(3),'clim',rng);
set(h(1),'ylim',rng);

end
