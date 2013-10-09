%Paprastuju iteraciju metoda naudosime ir lygciu sistemu sprendimui
% 8 savaite atsiskaitymas 

% 5 Variantas
% Daugianaris - x^3 - x^2 - 30x
% Metodai: skenavimo, Niutono
%
% funkcija g(x) = e^-pi * cos(x) * sin(x^2-1)   7 <= x <= 8
% Metodai: skenavimo, stygø

function L1
clc,close all
syms g x
%------------------------   PRADINIAI DUOMENYS  ---------------------------
% funkcijos simboline israiska
f='exp(-x) .* cos(x) .* sin(x.^2-1)'


range=[7,8] % TODO: parenkame saknu atskyrimo intervala 
grange = [-30, 30]


% eps=1e-9;  % TODO: parenkame tikslumo reiksme 
nitmax=100;% TODO: parenkame didziausia leistina iteraciju skaiciu
step = 0.001

% method='chords';


% braizomas funkcijos grafikas
npoints=1000; x=range(1): (range(2)-range(1))/(npoints-1) :range(2);  
figure(1); grid on; hold on;
str=[f,'=0;   Stygø metodas']; title(str);
plot(x,eval(f),'r-');
plot(range,[0 0],'b-');

%------------------------   SPRENDIMAS  -----------------------------------
%xn ir xn1 - intervalo galai
%prec - tikslumas
xn=range(1);xn1=range(2);prec=1;

xcur = xn;
x0 = -29;

sprendiniai = [];
intervalai = [ ; ];

pr = range(1);

while pr < range(2)
    sp = skenavimas(f, pr, eps, step, range(2));
    if ~isempty(sp) 
        sprendiniai = [sprendiniai sp];
    else
        break
    end
    int = [pr; sp + step * 10];
    disp('Rastas sprendinys intervale')
    disp(int)
    disp('\n')
    intervalai = [intervalai int];
    pr = sp + step;
end

isize = size(intervalai)

for i=1:isize(2)
    styguMetodas(f, intervalai(1, i) ,intervalai(2, i) ,prec,nitmax, sprendiniai)
end

plot(x,eval(f),'r-');plot(range,[0 0],'b-');plot(x0,0,'mp');h = findobj(gca,'Type','line');h1=h(1);

disp('.... Patikriname saknies reiksme, naudodami MATLAB funkcija fsolve: ....')
fprintf('\n')

    % fsolve iejimo pirmu parametru gali buti simboliu masyvas(char), kuriame irasyta lygties funkcijos israiska. 
    % Lygties funkcijos argumentas yra simbolis x. 
    % Antras parametras yra x reiksme. x gali buti vektorius.
    
 x0=-10; sprendinys = fsolve(char(f),x0); fprintf(1,'\npradinis artinys  x0=%g,  sprendinys=%g',x0,sprendinys);
 x0=5;   sprendinys = fsolve(char(f),x0); fprintf(1,'\npradinis artinys  x0=%g,  sprendinys=%g',x0,sprendinys);
 fprintf(1,'\n');
end

function spr = skenavimas(f,a, eps, step, max)
    prec = 1;
    b = a + step; 
    while (prec > eps) && (a ~= b) && (b < max)
        x = a; 
        fa = eval(f);
        x = b;
        fb = eval(f);  
        
       x = a;
        if sign(fa) == sign(fb) ~= 0
            a = b;
            b = b + step;
         else
             a = a + step / 10000;
        end 
 
        fx = eval(f);
        prec = abs(fx);
    end
    if prec <= eps
    spr = a;
    plot(spr, 0, 'mp');
    else spr = [];
    end
end

function styguMetodas(f, xn, xn1, prec, nitmax, sprendiniai)
nit=0;

while prec > eps
    %TODO kas yra h? Kas yra nit?
    
    %TODO bbz
    nit=nit+1;
    if nit > nitmax, fprintf('Virsytas leistinas iteraciju skaicius');break;end
    
    % Nupaisom intervalo galus
    plot(xn,0,'mp');h = findobj(gca,'Type','line');h1=h(1); % paskutinio grafinio objekto valdiklis irasomas handle masyvo priekyje
    plot(xn1,0,'cp');h = findobj(gca,'Type','line');h2=h(1);
     
            % randam f-jos reiksmes taskuose xn ir xn1
            x=xn;fxn=eval(f);
            x=xn1;
            fxn1=eval(f);
    
            % Apskaiciuojam koeficienta k (zr skaidres)
            k=abs(fxn/fxn1);
    
            % Apskaiciuojam xmid (zr skaidres)
            xmid=(xn+k*xn1)/(1+k);
            
            % Pazymim xmid grafike
            plot(xmid,0,'gs');
            
            plot([xn,xn1],[fxn,fxn1],'g-');
            h = findobj(gca,'Type','line');h3=h(1:2);
    
    % Apskaiciuojame funkcijos reiksme xmid taske
    x=xmid;
    fxmid=eval(f);
    
    % jeigu pradzioje tikriname kairi taska
    % Tikriname kairi intervalo taska
    x=xn;
    fxn=eval(f);
    
    % (zr skaidres)
    if sign(fxmid) == sign(fxn), xn=xmid;
    else xn1=xmid;
    end
    
%     x=xn1;fxn1=eval(f);
%     if sign(fxmid) == sign(fxn1), xn1=xmid;
%     else, xn=xmid;
%     end
        
    pause(1)
    delete(h1);delete(h2);delete(h3);
    
    x = xmid;
    prec=abs(fxmid); 
    % Idedam x i sprendiniu sarasa ir pazymim grafike
    if (prec <= eps && ~ismember(x, sprendiniai))
        %sprendiniai = [sprendiniai x]
        plot(x, 0, 'gs')
    end
    fprintf(1,'iteracija %d    tikslumas= %g \n',nit,prec);
    end
end