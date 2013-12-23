function plotmesh2d(edges,verts,ports,plot_normals)

vx = verts(:,1);
vy = verts(:,2);
plot(vx(edges)', vy(edges)', 'b-');
xlabel('x');
ylabel('y');

hold on
ec = edges(cell2mat(ports'),:);
plot(vx(ec)', vy(ec)', 'r-');
%plot(vx(ec)', vy(ec)', 'g-');

if plot_normals,
	% Normals
	
	r1 = verts(edges(:,1),:);
	r2 = verts(edges(:,2),:);
	c = (r1 + r2)/2;
	edges = r2 - r1;
	l = sqrt(sum(edges.^2,2));
	t = edges ./ l(:,ones(1,2));
	n = [ t(:,2) -t(:,1) ];
	
	quiver(c(:,1),c(:,2),n(:,1),n(:,2),0.3);
end

hold off
