%example = example1;
%example = example2;
%example = example3;
example = example4;









pde = example.pde;
option = example.option;
mesh = example.mesh;
h = mesh.h;
tau = mesh.tau;
T = mesh.T;
IT = mesh.IT;

Output(pde,option,IT,T,h,tau);















function Output(pde,option,IT,T,h,tau)
%% Show the result of solving the hyperbolic or conservation equation at time T in the interval IT[a,b]
% with the space step h and time step tau
% pde determine the exact solution to the equation, or the flux f or the
% convection term a (a is a constant or a=f'(u)) and the boundary with Fixed  or  'Periodic'
% option determine the solving method Upwind, LaxFried or  LaxWendroff and
% some contral of the output
Io = InitialTnertval(IT,pde.a,T,pde.bdtype);

xo = Io(1):h:Io(2);
idx = (xo<=IT(2)+1e-13)&(xo>=IT(1)-1e-13);
x = xo(idx);

Uo = pde.uo(xo);

O = size(option);

U = cell(O(1),2);

for k = 1:O(1)
    [Us,Cor] = Solve(Uo,T,h,tau,pde.a,option{k},pde.bdtype);
    U{k,1} = Us(idx); 
    U{k,2} = Cor;
end



if isfield(pde,'exactu')
    
    exactU = pde.exactu(x,T);
    Xlegend = cell(1,O(1));
    if isfield(pde,'ordershow')
        Os = O(1)+1;
    else
        Os = O(1);
    end
    for k = 1:O(1)
        if isfield(pde,'figureshow')
           subplot(Os,1,k);
            plot(x,U{k},'o');hold on
            plot(x,exactU,'-');
            legend(option{k},'exactu');
            title(['time =',num2str(T),', Cournum = ',num2str(U{k,2})]);
        end
        if isfield(pde,'ordershow')
           [rate,err,hs]=Rate(pde,option{k},T,h,tau,Io,pde.bdtype,pde.lp);
           subplot(Os,1,Os);
           plot(-log10(hs),log10(err),'-');hold on
           Xlegend{k}  = ['time = ',num2str(T),' ',option{k},', order = ',num2str(rate(end))];
           fprintf(['The err and order of ',option{k},' at time ',num2str(T),'\n']);
           fprintf('&h:\t')
           for fpi = 1:length(hs)
               fprintf('&%0.2e\t',hs(fpi));
           end
           fprintf('\\\\ \n');
           fprintf('&&err\t');
           for fpi = 1:length(err)
               fprintf('&%0.2e\t',err(fpi));
           end
           fprintf('\\\\ \n');
           fprintf('&&rate\t&-\t\t');
           for fpi = 1:length(rate)
               fprintf('&%0.2f\t\t',rate(fpi));
           end
           fprintf('\\\\ \n');
       end
    end
    if isfield(pde,'ordershow')
     subplot(Os,1,Os);
     legend(Xlegend);
     xlabel('-log10(h)');
     ylabel('log10(err)');
     title('order');
    end
else
     Xlegend = cell(1,O(1));
    for k = 1:O(1)
        plot(x,U{k});hold on
         Xlegend{k}  = [option{k},', Cournum = ',num2str(U{k,2})];
    end
   legend(Xlegend);
   title(['time = ',num2str(T)]);
end






end









function Io = InitialTnertval(I,a,T,bdtype)
%%determin the intial intervar 
%% a.Iu is a  priori estimate

switch bdtype %choice boundary type
    case 'Fixed'
        if isreal(a)
            Io = [I(1)-abs(a)*T-1,I(2)+abs(a)*T+1];
        else
            A = max(abs(a.a(a.Iu(1):(a.Iu(2)-a.Iu(1))/100:a.Iu(2))));
            Io = [I(1)-A*T-5,I(2)+A*T+5];
        end
    case 'Periodic'
        Io = I;
end

end



function [Us,Cor] = Solve(Uo,T,h,tau,a,option,bdtype)
switch option
    case 'upwind'
         [Us,Cor]= UpWindCIR(Uo,h,tau,a,T,bdtype);
    case 'LF'
         [Us,Cor]=  LaxFriedrichs(Uo,h,tau,a,T,bdtype);
    case 'LW'
         [Us,Cor] = LaxWendroff(Uo,h,tau,a,T,bdtype);
end
end
       


        

    