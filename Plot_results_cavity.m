% Visualization of modeling results of 2D cavity flows.
% Cavity model was developed with Fortran + MPI, while this script 
% just plot the results


clear 
tic

result = load('final_results_mpi.txt');

x = result(:,1);
y = result(:,2);
p = result(:,3);
u = result(:,4);
v = result(:,5);


[X,Y] = meshgrid(unique(x),unique(y));
%{
p_map=griddata(x,y,p,X,Y);
u_map=griddata(x,y,u,X,Y);
v_map=griddata(x,y,v,X,Y);
p_map_ext=griddata(x,y,pexact,X,Y);
u_map_ext=griddata(x,y,uexact,X,Y);
v_map_ext=griddata(x,y,vexact,X,Y);
%}
%{
[~,p_map]=meshgrid(x,p);
[~,u_map]=meshgrid(x,u);
[~,v_map]=meshgrid(x,v);
[~,p_map_ext]=meshgrid(x,pexact);
[~,u_map_ext]=meshgrid(x,uexact);
[~,v_map_ext]=meshgrid(x,vexact);
%}
p_map = reshape(p,[length(unique(x)),length(unique(y))])';%so that x as row, y as column
u_map = reshape(u,[length(unique(x)),length(unique(y))])';
v_map = reshape(v,[length(unique(x)),length(unique(y))])';

toc
%%
tic
figure
contourf(X,Y,p_map)
colorbar
axis([0,max(x),0,max(y)])
set(gca,'fontsize',15)
axis equal
title('p ','fontsize',14)
%axis equal

figure
contourf(X,Y,u_map)
colorbar
axis([0,max(x),0,max(y)])
set(gca,'fontsize',15)
axis equal
title('u (m/s)','fontsize',14)
%axis equal

figure
contourf(X,Y,v_map)
colorbar
axis([0,max(x),0,max(y)])
set(gca,'fontsize',15)
axis equal
title('v (m/s)','fontsize',14)
%axis equal
%%
scrsz = get(0,'ScreenSize');
field_u = zeros( size(u_map) );
field_v = zeros( size(v_map) );
for i = 1 : sqrt( size(u,1) )
    for j = 1 : sqrt(size(u,1))
        field_u(i,j) = u_map(i,j);
        field_v(i,j) = v_map(i,j);
    end
end

amp = sqrt( field_u.^2+field_v.^2 );
field_u = field_u./amp;
field_v = field_v./amp;

figure
hold on
contourf(X,Y,u_map)
axis equal
colorbar
quiver(X,Y,field_u,field_v,1,'k')
axis([0,max(x),0,max(y)])
set(gca,'fontsize',15)
title('Velocity Field','fontsize',18)

set(gcf,'Position',[1 scrsz(4) scrsz(3)*0.5 scrsz(4)*0.7])
set(gcf,'units',get(gcf,'paperunits'));
set(gcf,'paperposition',get(gcf,'position'));
saveas(gcf,'Velocity_Field_Re100.png')

toc
