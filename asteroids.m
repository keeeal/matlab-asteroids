function [ ] = asteroids(start_level, varargin)
% ASTEROIDS is a matlab version of the classic game, asteroids. Use the arrow
%   keys to navigate your ship and spacebar to shoot. Destroy as many
%   asteroids as you can without being hit by them.
%
%   James Keal
%   21/9/13

close all
if nargin < 1
    start_level = 1;
end

% initialise figure
scrsz = get(0,'screensize');
wdth = scrsz(4)/2;
winsz = [scrsz(3)/2-wdth/2,wdth/2,wdth,wdth];
figure('units','pixels',...
    'name','Asteroids',...
    'menubar','none',...
    'numbertitle','off',...
    'position',winsz,...
    'color','k',...
    'keypressfcn',@kpf,...
    'closereq',@crf,...
    'busyaction','cancel',...
    'renderer','opengl');

% initialise ship variables
ship_x = [-1,1,0,-1];
ship_y = [-1,-1,2,-1];
s_pos = [0,0];
s_vel = [0,0];
s_rot = 0;

% initialise asteroid variables
a = 0;
a_pos_x = []; a_pos_y = [];
ast_x = zeros(13,1);
ast_y = zeros(13,1);
a_vel_x = []; a_vel_y = [];
a_rad = [];
a_delay = 100;
a_col = 'w';

% initialise bullet variables
b = 0;
bul_x = []; bul_y = [];
b_vel_x = []; b_vel_y = [];
b_time = [];
b_col = 'c.';

% initialise laser variables
las_x = [0,0]; las_y = [0,0];
l_time = 0;
l_col = 'r';

% initialise powerup variables
pu_delay = 400;
pu_x = []; pu_y = [];
pu_arr = [0,0,0,0,0];
pu_lim = [0,0,0,0,0];
pu_clr = ['y','g','b','m','r'];
nxt_pu = 1;
pu_str = 'pea-shooter';

% initialise other variables
score = 0;
cheat = zeros(7,1);
s_col = 'w';
level = start_level;
level_str = int2str(level);
endgame = 0;

% generates a random-shaped asteroid
    function [] = asteroid(x,y,vx,vy,r)
        a = a + 1;
        a_rad(a) = r;
        a_pos_x(a) = x;
        a_pos_y(a) = y;
        ast_x(:,a) = [(r+r*rand(12,1)).*sin((0:11)'*pi/6);0]+x;
        ast_x(end,a) = ast_x(1,a);
        ast_y(:,a) = [(r+r*rand(12,1)).*cos((0:11)'*pi/6);0]+y;
        ast_y(end,a) = ast_y(1,a);
        a_vel_x(a) = vx;
        a_vel_y(a) = vy;
    end

% destroy an asteroid
    function [] = split(logic_a)
        del_r = a_rad(~logic_a);
        del_px = a_pos_x(~logic_a);
        del_py = a_pos_y(~logic_a);
        a_pos_x = a_pos_x(logic_a);
        a_pos_y = a_pos_y(logic_a);
        a_vel_x = a_vel_x(logic_a);
        a_vel_y = a_vel_y(logic_a);
        a_rad = a_rad(logic_a);
        logic_a = logical(logic_a.*(1:a));
        ast_x = ast_x(:,logic_a);
        ast_y = ast_y(:,logic_a);
        a = length(a_pos_x);
        if sum(del_r)
            score = score + floor(4000/sum(del_r));
        end
        
        % split asteroid
        if del_r > 2
            nv = 2*rand(1,2)-1;
            asteroid(del_px(1),del_py(1),nv(1),-nv(2),del_r(1)/2);
            asteroid(del_px(1),del_py(1),-nv(1),nv(2),del_r(1)/2);
        end
    end

% calculates edge-wrap correction
    function [w_out] = wrap(w_in)
        w_out = 200*sum([(w_in<-100);-(w_in>100)]);
    end

% increase the level
    function [] = levelup
        if level < 12
            level = level + 1;
            level_str = int2str(level);
        end
    end

% set powerup
    function [] = setpu(pu_num)
        if pu_num == 0
            pu_arr = pu_arr-pu_arr;
            pu_lim = [0,0,0,3,0];
            b_col = 'c.';
            pu_str = 'pea-shooter';
        elseif pu_num <= length(pu_arr)
            pu_arr(pu_num) = 1;
            switch pu_num
                case 1
                    b_col = 'yd';
                    pu_str = 'shotgun';
                case 2
                    b_col = 'gx';
                    pu_str = 'tri-shield';
                case 3
                    b_col = 'b*';
                    pu_str = 'long-shot';
                case 4
                    b_col = 'mo';
                    pu_str = 'radiation';
                case 5
                    pu_str = 'laser';
            end
        end
    end

% start timer
tmr = timer('Name','the_timer',...
    'Period',0.06,...
    'StartDelay',1,...
    'ExecutionMode','fixedrate',...
    'TimerFcn',@replot);
start(tmr)

% draws updated plot
    function [] = replot(varargin)
        % refresh ship
        s_wrap = wrap(s_pos);
        ship_x = ship_x + s_vel(1) + s_wrap(1);
        ship_y = ship_y + s_vel(2) + s_wrap(2);
        s_pos = s_pos + s_vel + s_wrap;
        
        % refresh bullets
        bt = logical(b_time);
        b_vel_x = b_vel_x(bt);
        b_vel_y = b_vel_y(bt);
        bul_x = bul_x(bt) + b_vel_x + wrap(bul_x(bt));
        bul_y = bul_y(bt) + b_vel_y + wrap(bul_y(bt));
        b_time = b_time(bt) - 1;
        b = length(b_time);
        
        % refresh laser
        if l_time
            l_time = l_time - 1;
        else
            las_x = [0,0];
            las_y = [0,0];
        end
        
        % power ups with autofire
        if pu_arr(2)
               for iii = 1:3
               shoot(s_rot+(2*iii+1)*pi/3,8);
               end
        end
        if pu_arr(4) && pu_lim(4)
            for iii = 0:pi/64:2*pi
                shoot(iii,150);
            end
            pu_lim(4) = pu_lim(4)-1;
        end
        
        % refresh asteroids
        a_delay = a_delay - 1;
        if a
            a_wrap_x = (1+a_rad/50).*wrap(a_pos_x-2*sign(a_pos_x).*a_rad);
            a_wrap_y = (1+a_rad/50).*wrap(a_pos_y-2*sign(a_pos_y).*a_rad);
            [xj,~] = meshgrid(a_vel_x+a_wrap_x,(0:12));
            ast_x = ast_x + xj;
            [yj,~] = meshgrid(a_vel_y+a_wrap_y,(0:12));
            ast_y = ast_y + yj;
            a_pos_x = a_pos_x + a_vel_x + a_wrap_x;
            a_pos_y = a_pos_y + a_vel_y + a_wrap_y;
            
            % collisions
            if b
                out_bul = zeros(b,a);
                for j = 1:a
                    out_bul(:,j) = ~inpolygon(bul_x,bul_y,ast_x(:,j),ast_y(:,j));
                end
                
                % remove bullet
                prod_b = logical(prod(out_bul,2));
                b_vel_x = b_vel_x(prod_b);
                b_vel_y = b_vel_y(prod_b);
                bul_x = bul_x(prod_b);
                bul_y = bul_y(prod_b);
                b_time = b_time(prod_b);
                
                % remove asteroid
                prod_a = logical(prod(out_bul,1));
                split(prod_a);
            end
            
            if l_time
                out_las = zeros(1,a);
                for j = 1:a
                    out_las(j) = sum(~~polyxpoly(las_x,las_y,ast_x(:,j),ast_y(:,j)));
                end
                
                % remove asteroid
                prod_a = logical(~out_las);
                split(prod_a);
            end
        end
        
        % generate asteroids
        if a_delay < 0
            in_r = 10*rand + 5;
            in_x = rand;
            if in_x < 0.25
                in_y = (8*in_x-1)*(100+2*in_r);
                in_x = 100+2*in_r;
            elseif in_x < 0.5
                in_y = (8*in_x-3)*(100+2*in_r);
                in_x = -100-2*in_r;
            elseif in_x < 0.75
                in_y = 100+2*in_r;
                in_x = (8*in_x-5)*(100+2*in_r);
            else
                in_y = -100-2*in_r;
                in_x = (8*in_x-7)*(100+2*in_r);
            end
            asteroid(in_x,in_y,-in_x/100,-in_y/100,in_r);
            a_delay = floor(rand*(250-20*level));
        end
        
        % generate powerup square
        pu_s = sign(pu_delay);
        pu_delay = pu_delay - pu_s;
        if ~pu_delay
            if pu_s > 0
                pu_pos = 200*rand(1,2)-100;
                pu_x = [pu_pos(1)-4,pu_pos(1)+4,pu_pos(1)+4,pu_pos(1)-4,pu_pos(1)-4];
                pu_y = [pu_pos(2)-4,pu_pos(2)-4,pu_pos(2)+4,pu_pos(2)+4,pu_pos(2)-4];
                pu_delay = -200;
                nxt_pu = floor(rand*length(pu_clr)+1);
                setpu(0);
            else
                pu_x = []; pu_y = [];
                pu_delay = 400;
                levelup;
            end
        end
        
        if pu_delay < 0
            if inpolygon(ship_x,ship_y,pu_x,pu_y)
                setpu(nxt_pu);
            end
        end
        
        % plot
        plot(ship_x,ship_y,s_col,ast_x,ast_y,a_col,bul_x,bul_y,b_col,las_x,las_y,l_col,pu_x,pu_y,pu_clr(nxt_pu));
        axis([-100,100,-100,100]);
        set(gca,'position',[0,0,1,1],'visible','off');
        score = score + 17;
        for i = 1:a
            if inpolygon(s_pos(1),s_pos(2),ast_x(:,i),ast_y(:,i))
                if cheat(7)
                    for k = 0:pi/32:2*pi
                        shoot(k);
                    end
                else
                    text(-66,10,'GAME OVER','Color','w','FontSize',50);
                    text(-40,-10,['          score: ',int2str(score),char(10),'press return to play again'],'Color','w','FontSize',15);
                    file = 'scores_(asteroids).m';
                    if exist(file, 'file') == 2
                        load(file,'-mat','scores');
                    else
                        scores = (1e6:-1e5:1e5)';
                    end
                    scores = sort([scores;score],'descend');
                    scores = scores(1:end-1);
                    save(file,'-mat','scores');
                    your_score = score
                    high_scores = scores
                    endgame = 1;
                    stop(tmr);
                    return
                end
            end
        end
        text(80,-95,int2str(score),'Color','w');
        text(-90,-95,['level ',level_str],'Color','w');
        text(-90,95,pu_str,'Color','w');
    end

% key press function
    function [] = kpf(varargin)
        if strcmp(varargin{2}.Key,'uparrow')
            s_vel = s_vel + 0.25*[sin(s_rot),cos(s_rot)];
        elseif strcmp(varargin{2}.Key,'rightarrow')
            turn(0.4);
        elseif strcmp(varargin{2}.Key,'leftarrow')
            turn(-0.4);
        elseif strcmp(varargin{2}.Key,'downarrow')
            reverse;
        elseif strcmp(varargin{2}.Key,'space')
            if pu_arr(1)
                shoot(s_rot-0.15);
                shoot(s_rot);
                shoot(s_rot+0.15);
            elseif pu_arr(3)
                shoot(s_rot,200);
            elseif pu_arr(5)
                laser(s_rot);
            else
                shoot(s_rot);
            end
            
            % godmode
        elseif strcmp(varargin{2}.Key,'g')
            cheat(1) = 1;
        elseif strcmp(varargin{2}.Key,'o') && prod(cheat(1:4))
            cheat(5) = 1;
        elseif strcmp(varargin{2}.Key,'d') && prod(cheat(1:5))
            cheat(6) = 1;
        elseif strcmp(varargin{2}.Key,'m') && prod(cheat(1:3))
            cheat(4) = 1;
        elseif strcmp(varargin{2}.Key,'o') && cheat(1)
            cheat(2) = 1;
        elseif strcmp(varargin{2}.Key,'d') && prod(cheat(1:2))
            cheat(3) = 1;
        elseif strcmp(varargin{2}.Key,'e') && prod(cheat(1:6))
            cheat(7) = 1;
            s_col = 'r';
        end
        
        if cheat(7)
            if strcmp(varargin{2}.Key,'q')
                levelup;
            elseif strcmp(varargin{2}.Key,'a')
                if level > 1
                    level = level-1;
                    level_str = int2str(level);
                end
            else
                num_in = str2double(varargin{2}.Key);
                if ~isempty(num_in)
                    setpu(num_in);
                end
            end
        end
        
        if endgame
            if strcmp(varargin{2}.Key,'return')
                asteroids(start_level);
            end
        end
    end

    function [] = turn(dir)
        s_rot = mod(s_rot+dir,2*pi);
        x1 = sqrt(2)*sin(-3*pi/4+s_rot);
        ship_x = [x1,sqrt(2)*sin(3*pi/4+s_rot),2*sin(s_rot),x1]+s_pos(1);
        y1 = sqrt(2)*cos(-3*pi/4+s_rot);
        ship_y = [y1,sqrt(2)*cos(3*pi/4+s_rot),2*cos(s_rot),y1]+s_pos(2);
    end

    function [] = reverse
        rev_dir =  pi+sign(s_vel(1))*acos(s_vel(2)/norm(s_vel));
        rev = rev_dir-s_rot;
        if abs(rev) > 0.4
            if rev > s_rot
                turn(0.4);
            else
                turn(-0.4);
            end
        elseif rev
            turn(rev);
        end
    end

    function [] = shoot(dir,time)
        b = b + 1;
        bul_x(b) = s_pos(1);
        bul_y(b) = s_pos(2);
        b_vel_x(b) = 3*sin(dir);
        b_vel_y(b) = 3*cos(dir);
        b_time(b) = 25;
        if nargin == 2
            b_time(b) = time;
        end
    end

    function [] = laser(dir)
        las_x = [s_pos(1),200*sin(dir)];
        las_y = [s_pos(2),200*cos(dir)];
        l_time = 5;
    end

% close request funtion
    function [] = crf(varargin)
        stop(tmr); delete(tmr);
        delete(varargin{1});
    end

end
