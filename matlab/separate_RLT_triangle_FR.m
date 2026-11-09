function [Bcell, l, viol_vect] = separate_RLT_triangle_FR(Y, eps, max_sep_ineq, max_new_ineq)

    n = size(Y, 1);
    
    Bcell = cell(1, max_sep_ineq);
    l = zeros(max_sep_ineq, 1);
    viol_vect = zeros(max_sep_ineq, 1);
    c = 1;
    
    stop = false;
    for i=1:n
        for j=i+1:n
            % Yii-Yij >= 0
            viol = Y(i,i)-Y(i,j);
            if viol <= -eps
                Bcell{c} = sparse([i, i, j], [i, j, i], [1.0, -0.5, -0.5], n, n);
                l(c) = 0;
                viol_vect(c) = viol;
                c = c + 1;
            end
            % Yij-Yii-Yjj >= -1
            viol = Y(i,j)-Y(i,i)-Y(j,j)+1;
            if viol <= -eps
                Bcell{c} = sparse([i, j, i, j], [j, i, i, j], [0.5, 0.5, -1.0, -1.0], n, n);
                l(c) = -1;
                viol_vect(c) = viol;
                c = c + 1;
            end
            % Yjj-Yij >= 0
            viol = -Y(i,j)+Y(j,j);
            if viol <= -eps
                Bcell{c} = sparse([j, i, j], [j, j, i], [1.0, -0.5, -0.5], n, n);
                l(c) = 0;
                viol_vect(c) = viol;
                c = c + 1;
            end
            for t=j+1:n
                 if (i ~= t)
                     % -Yii-Yjj-Ytt+Yij+Yit+Yjt >= -1
                     viol = -Y(i,i)-Y(j,j)-Y(t,t)+Y(i,j)+Y(i,t)+Y(j,t) + 1;
                     if viol <= -eps
                         Bcell{c} = sparse( ...
                             [i, j, t, i, j, i, t, j, t], ...
                             [i, j, t, j, i, t, i, t, j], ...
                             [-1.0, -1.0, -1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], n, n);
                         l(c) = -1;
                         viol_vect(c) = viol;
                         c = c + 1;
                     end
                     % -Yij-Yit+Yii+Yjt >= 0
                     viol = -Y(i,j)-Y(i,t)+Y(i,i)+Y(j,t);
                     if viol <= -eps
                         Bcell{c} = sparse( ...
                             [i, j, i, t, i, j, t], ...
                             [j, i, t, i, i, t, j], ...
                             [-0.5, -0.5, -0.5, -0.5, 1.0, 0.5, 0.5], n, n);
                         l(c) = 0;
                         viol_vect(c) = viol;
                         c = c + 1;
                     end
                     % -Yij-Yjt+Yjj+Yit >= 0
                     viol = -Y(i,j)-Y(j,t)+Y(j,j)+Y(i,t);
                     if viol <= -eps
                         Bcell{c} = sparse( ...
                             [i, j, j, t, j, i, t], ...
                             [j, i, t, j, j, t, i], ...
                             [-0.5, -0.5, -0.5, -0.5, 1.0, 0.5, 0.5], n, n);
                         l(c) = 0;
                         viol_vect(c) = viol;
                         c = c + 1;
                     end
                     % -Yit-Yjt+Ytt+Yij >= 0
                     viol = -Y(i,t)-Y(j,t)+Y(t,t)+Y(i,j);
                     if viol <= -eps
                         Bcell{c} = sparse( ...
                             [i, t, j, t, t, i, j], ...
                             [t, i, t, j, t, j, i], ...
                             [-0.5, -0.5, -0.5, -0.5, 1.0, 0.5, 0.5], n, n);
                         l(c) = 0;
                         viol_vect(c) = viol;
                         c = c + 1;
                     end
                     
                     if c >= max_sep_ineq
                         stop=true;
                         break;
                     end
                     
                 end
            end
            
            if stop
                break;
            end
            
        end
        
        if stop
            break;
        end
        
    end
    
    n_ineq = c - 1;
    if n_ineq <= max_new_ineq
        Bcell = Bcell(1:n_ineq);
        l = l(1:n_ineq);
    else    
        viol_vect = abs(viol_vect(1:n_ineq));
        [~, id_sorted] = sort(viol_vect, 'descend');
        Bcell = Bcell(id_sorted);
        l = l(id_sorted);
        Bcell = Bcell(1:max_new_ineq);
        l = l(1:max_new_ineq);
    end
    
end