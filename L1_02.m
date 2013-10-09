function L1_02
syms f x

f = x^3 - x^2 - 30 * x
df = diff(f,x)
range = [-7.5, 7.5]
x0 = -25;
nitmax = 100;
beta = 1;

eps = 1e-9;
pr = range(1);
step = 0.1

figure(1); grid on; hold on; axis equal; str=[char(f),'=0;   Niutono metodas'];title(str);
npoints=1000;
xrange = range(1): (range(2)-range(1))/(npoints-1) :range(2);
x = xrange;

plot(x,eval(f),'r-');
plot(range,[0 0],'b-');
plot(x0,0,'mp');
h = findobj(gca,'Type','line');
h1=h(1);

intervalai = [];
sprendiniai = [];
while pr < range(2)
    sp = skenavimas(f, pr, 1e-3, step, range(2));
    if ~isempty(sp) 
        sprendiniai = [sprendiniai sp];
    else
        break
    end
    int = [pr; sp + 1];
    disp('Rastas sprendinys intervale')
    disp(int)
    disp('\n')
    intervalai = [intervalai int];
    pr = sp + step;
end

isize = size(intervalai)

for i=1:isize(2)
    Newton(f, df, intervalai(1,i) + 1);
end

    Newton(f, df, 5);

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
             a = a + step / 100;
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

function Newton(f, df, x0)
xn=x0;prec=1;nit=0;
while prec > eps    % iteracijos 
    nit=nit+1;
    
    if nit > nitmax, fprintf('Virsytas leistinas iteraciju skaicius');
        break;
    end
    
            x=xn;
            fxn=eval(f);
            dfxn=eval(df);
            xn1=xn-beta*fxn/dfxn;   % daugikliu alpha galima apriboti x prieaugi, tikintis
                                     % geresnio konvergavimo   
            plot([xn,xn,xn1],[0,fxn,0],'g-');
            delete(h1);plot(xn1,0,'mp');h = findobj(gca,'Type','line');h1=h(1);
            xn=xn1; 
    x=xn;
    fxn=eval(f);
    prec=abs(fxn);
    fprintf(1,'iteracija %d  x= %g  prec= %g \n',nit,xn,prec);
end
plot(xn,fxn,'k*');plot(xn,fxn,'ko');
xn
nit

end

end