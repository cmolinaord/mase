function [V_x_e,V_x_w,V_y_n,V_y_s] = SH_vel(Xpos,Ypos,w)

	V_x_e = zeros(w.ny,w.nx); % X axis velocity (east face)
	V_x_w = zeros(w.ny,w.nx); % X axis velocity (west face)
	V_y_n = zeros(w.ny,w.nx); % Y axis velocity (north face)
	V_y_s = zeros(w.ny,w.nx); % Y axis velocity (south face)

	V_x = zeros(w.ny,w.nx); % X axis velocity (node P)
	V_y = zeros(w.ny,w.nx); % Y axis velocity (node P)

	for i = 2:w.nx-1
		for j = 2:w.ny-1

			deltaXe = (Xpos(i+1)-Xpos(i));
			deltaXw = (Xpos(i)-Xpos(i-1));
			deltaYn = (Ypos(j+1)-Ypos(j));
			deltaYs = (Ypos(j)-Ypos(j-1));

			V_x_e(j,i) = 2*Ypos(j)*(1-(Xpos(i)+deltaXe/2)^2); % X axis velocity (east face)
			V_x_w(j,i) = 2*Ypos(j)*(1-(Xpos(i)-deltaXw/2)^2); % X axis velocity (west face)
			V_y_n(j,i) = -2*Xpos(i)*(1-(Ypos(j)+deltaYn/2)^2); % Y axis velocity (north face)
			V_y_s(j,i) = -2*Xpos(i)*(1-(Ypos(j)-deltaYs/2)^2); % Y axis velocity (south face)
			V_x(j,i) = 2*Ypos(j)*(1-Xpos(i)^2); % X axis velocity (node)
			V_y(j,i) = -2*Xpos(i)*(1-Ypos(j)^2); % Y axis velocity (node)
		end
	end
end
