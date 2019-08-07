function new_fig=inset(fig_num,main_handle, inset_handle,inset_size)

% The function plotting figure inside figure (main and inset) from 2 existing figures.
% inset_size is the fraction of inset-figure size, default value is 0.35
% The outputs are the axes-handles of both.
% 
% 
% Moshe Lindner, August 2010 (C).

if nargin==3
    inset_size=0.35;
end

inset_size=inset_size*.7;
figure(fig_num)
new_fig=gcf;
main_fig = findobj(main_handle,'Type','axes');
h_main = copyobj(main_fig,new_fig);
set(new_fig,'Position',get(main_handle,'Position'))
inset_fig = findobj(inset_handle,'Type','axes');
h_inset = copyobj(inset_fig,new_fig);
ax=get(main_fig,'Position');
set(h_inset,'Position', [ax(1)+ax(3)-inset_size/2 .5*ax(2)+ax(4)-inset_size inset_size/2 inset_size])
close(main_handle)
close(inset_handle)