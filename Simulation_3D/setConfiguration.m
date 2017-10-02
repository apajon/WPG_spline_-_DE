function [move_dirichlet,first_opt,trajectory,material,w,s1,s2,s3,p_ini_v] = setConfiguration(type,l_spli)
    s1 = 'results/CMAME/';
    s3 = '/';
    s4 = 'polynome';
    s5 = '.mat';    
    if type == 1
        s2 = 'Neoprene';
        move_dirichlet = 1;
        first_opt = 1;
        trajectory = 1;
        material = 1;
        w = 1;
    elseif type == 2
        s2 = 'ButylRubber';
        move_dirichlet = 1;
        first_opt = 1;
        trajectory = 1;
        material = 2;
        w = 1;
    elseif type == 3
        s2 = 'Straight';
        move_dirichlet = 1;
        first_opt = 1;
        trajectory = 2;
        material = 1;
        w = 1;
    elseif type == 4       
        s2 = '1-1Neoprene';
        s = strcat(s1,s2,s3,s4,s2,s5);
        move_dirichlet = 1;
        first_opt = 0;
        trajectory = 1;
        material = 1;
        w = 1;
    elseif type == 5
        s2 = '1-1000Neoprene';
        s = strcat(s1,s2,s3,s4,s2,s5);
        move_dirichlet = 1;
        first_opt = 0;
        trajectory = 1;
        material = 1;
        w = 1000;
    elseif type == 6
        s2 = '1-1ButylRubber';
        s = strcat(s1,s2,s3,s4,s2,s5);   
        move_dirichlet = 1;
        first_opt = 0;
        trajectory = 1;
        material = 2;
        w = 1;
    elseif type == 7
        s2 = '1-1Straight';
        s = strcat(s1,s2,s3,s4,s2,s5); 
        move_dirichlet = 1;
        first_opt = 0;
        trajectory = 2;
        material = 1;
        w = 1;
    elseif type == 8
        s2 = '1-1NoDir';
        s = strcat(s1,s2,s3,s4,s2,s5);  
        move_dirichlet = 0;
        first_opt = 0;
        trajectory = 1;
        material = 1;
        w = 1;
    end
    if first_opt==1 
        if move_dirichlet==1
            p_ini_v = [ones(l_spli*l_spli,1); ...
                ones(l_spli-2,1); ones(l_spli-2,1);1;1];
        else
            p_ini_v = ones(l_spli*l_spli,1);
        end
    else
        a = load(s);
        p_ini_v = a.p_ini_v;
    end
end