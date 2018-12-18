function [Ke, KG, R_, Re_] = stiffness_matrix(data, Tmat, Tn, mat, dat, A, L, I_z, I_y, J)

    for e=1:data.n_el

        if Tmat(e,1) == 1
            E(e) = mat(1,1);
            G(e) = mat(1,2);
        elseif Tmat(e,1) == 2
            E(e) = mat(2,1);
            G(e) = mat(2,2);
        end

        alpha = dat(e,2);
        beta  = dat(e,3);
        gamma = dat(e,4);

        Ke_p = [E(e)*A(e)/L(e) 0 0 0 0 0 -E(e)*A(e)/L(e) 0 0 0 0 0;
                0 12*E(e)*I_z(e)/(L(e)^3) 0 0 0 6*E(e)*I_z(e)/(L(e)^2) 0 -12*E(e)*I_z(e)/(L(e)^3) 0 0 0 6*E(e)*I_z(e)/(L(e)^2);
                0 0 12*E(e)*I_y(e)/(L(e)^3) 0 -6*E(e)*I_y(e)/(L(e)^2) 0 0 0 -12*E(e)*I_y(e)/(L(e)^3) 0 -6*E(e)*I_y(e)/(L(e)^2) 0;
                0 0 0 G(e)*J(e)/L(e) 0 0 0 0 0 -G(e)*J(e)/L(e) 0 0;
                0 0 -6*E(e)*I_y(e)/(L(e)^2) 0 4*E(e)*I_y(e)/L(e) 0 0 0 6*E(e)*I_y(e)/(L(e)^2) 0 2*E(e)*I_y(e)/L(e) 0;
                0 6*E(e)*I_z(e)/(L(e)^2) 0 0 0 4*E(e)*I_z(e)/L(e) 0 -6*E(e)*I_z(e)/(L(e)^2) 0 0 0 2*E(e)*I_z(e)/L(e);
                -E(e)*A(e)/L(e) 0 0 0 0 0 E(e)*A(e)/L(e) 0 0 0 0 0;
                0 -12*E(e)*I_z(e)/L(e)^3 0 0 0 -6*E(e)*I_z(e)/L(e)^2 0 12*E(e)*I_z(e)/L(e)^3 0 0 0 -6*E(e)*I_z(e)/L(e)^2;
                0 0 -12*E(e)*I_y(e)/L(e)^3 0 6*E(e)*I_y(e)/L(e)^2 0 0 0 12*E(e)*I_y(e)/L(e)^3 0 6*E(e)*I_y(e)/L(e)^2 0;
                0 0 0 -G(e)*J(e)/L(e) 0 0 0 0 0 G(e)*J(e)/L(e) 0 0;
                0 0 -6*E(e)*I_y(e)/L(e)^2 0 2*E(e)*I_y(e)/L(e) 0 0 0 6*E(e)*I_y(e)/L(e)^2 0 4*E(e)*I_y(e)/L(e) 0;
                0 6*E(e)*I_z(e)/L(e)^2 0 0 0 2*E(e)*I_z(e)/L(e) 0 -6*E(e)*I_z(e)/L(e)^2 0 0 0 4*E(e)*I_z(e)/L(e)];

        R = [cos(beta)*cos(gamma), cos(beta)*sin(gamma), sin(beta);
            -sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma), -sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma), sin(alpha)*cos(beta);
            -cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma), -cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma), cos(alpha)*cos(beta)];

        Re = [R(1,1) R(1,2) R(1,3) 0 0 0 0 0 0 0 0 0;
              R(2,1) R(2,2) R(2,3) 0 0 0 0 0 0 0 0 0;
              R(3,1) R(3,2) R(3,3) 0 0 0 0 0 0 0 0 0;
              0 0 0 R(1,1) R(1,2) R(1,3) 0 0 0 0 0 0;
              0 0 0 R(2,1) R(2,2) R(2,3) 0 0 0 0 0 0;
              0 0 0 R(3,1) R(3,2) R(3,3) 0 0 0 0 0 0;
              0 0 0 0 0 0 R(1,1) R(1,2) R(1,3) 0 0 0;
              0 0 0 0 0 0 R(2,1) R(2,2) R(2,3) 0 0 0;
              0 0 0 0 0 0 R(3,1) R(3,2) R(3,3) 0 0 0;
              0 0 0 0 0 0 0 0 0 R(1,1) R(1,2) R(1,3);
              0 0 0 0 0 0 0 0 0 R(2,1) R(2,2) R(2,3);
              0 0 0 0 0 0 0 0 0 R(3,1) R(3,2) R(3,3)];

        K = Re' * Ke_p * Re;

        %store element matrix
        for r = 1:data.nnod * data.nd * 2
            for s = 1:data.nnod * data.nd * 2
                Ke(r,s,e) = K(r,s); % Stiffness matrix for each element
                Re_(r,s,e) = Re(r,s);
            end
        end

        for r = 1:3
            for s= 1:3
            	R_(r,s,e)  = R(r,s);
            end
        end

    end


    KG = zeros(data.ndof,data.ndof);

    for e = 1:data.n_el
        for i = 1:data.nnod * data.ngl
            I = Tn(e,i);
            for j = 1:data.nnod * data.ngl
                Je = Tn(e,j);
                KG(I,Je) = KG(I,Je) + Ke(i,j,e); % Global stiffness matrix
            end
        end
    end

end
