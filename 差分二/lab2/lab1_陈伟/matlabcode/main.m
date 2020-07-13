clear
global g_D f delta alpha Ftype;
alpha = 0.6;
pde = example3;
g_D = pde.g_D;
f = pde.f;
 delta = 10;
 delta_max = 10;
 %% 迭代方法
method = 'Adaptive_Newtow_delta';%'Netwon';%'Euler';'Adaptive_Newtow_delta'

%% F---正则或非正则化法
Ftype = 'MAWS_delta';%'MAWS'; %'MAWS_delta'

Nn = [16,32,64];
WideNn = [1,2,3,4];
% Nn = 32;
%  WideNn = 1;
%Test = 2;
Err = zeros(length(WideNn),length(Nn));%---------误差
K = zeros(length(WideNn),length(Nn));%---------迭代次数
Time =  zeros(length(WideNn),length(Nn));%---------消耗时间
Delta = zeros(length(WideNn),length(Nn));%--------终止时delta值

%% 初值


for in = 1:length(Nn)
    N = Nn(in);
    
    for iwiden = 1:length(WideNn)
            
            WideN = WideNn(iwiden);
            h = 1/N;
             [i,j] = GlobaltoLocalidx((1:(N-1)^2)',N);
             u = pde.u(i*h,j*h);


             G = G_theta(WideN);

            
            U = initial(i,j,h);
 %              U = rand(N-1);
 %            U = ones(N-1);
             U = U(:);
             k = 0;
             ks = 0;



             %% 
             tic
             switch method
                 case 'Euler'

                      Us = Eulerexplicit(U,h,G);

                     err = norm(U-Us,inf);
                     while and(k<5000,err>1e-8)
                             k = k+1;
                             U = Us;
                            Us = Eulerexplicit(U,h,G);

                             err = norm(U-Us,inf);
                             err
                     end

                 case 'Netwon'
                     Us = Newton_method(U,h,G);
                     err = norm(U-Us,inf)
                      while and(k<1000,err>h^4)
                             k = k+1;
                             U = Us;
                            Us = Newton_method(U,h,G);

                             err = norm(U-Us,inf);
                             err;
                      end

                 case 'Adaptive_Newtow_delta'
                     Ftype = 'MAWS_delta';
                     delta = delta_max;
                     Uo = U;
                     U = Newton_method(Uo,h,G);
                     errs = norm(Uo-U,inf);
                     while and(delta>1e-6,and(ks<N,errs>1e-8))


                         k_o = 0;
                         ks = ks + 1;

                         Us = Newton_method(U,h,G);
                         err = norm(U-Us,inf);
                          while and(k_o<1000,err>1e-8)
                                 k = k+1;
                                 k_o = k_o+1;
                                 U = Us;
                                 Us = Newton_method(U,h,G);

                                 err = norm(U-Us,inf);
                                 err
                          end
                          delta = delta/2;
                          errs = norm(Uo-U,inf)
                          Uo = U;
                     end  
             end
             K(iwiden,in) =  K(iwiden,in) + k;
             Err(iwiden,in) = Err(iwiden,in) + norm(u-U,inf);
             Delta(iwiden,in) = delta;
             Time(iwiden,in) = Time(iwiden,in) + toc;
    end
end
 %norm(u-U,inf)
 
 
 %% 打印结果
 
 for i = 1:size(K,1)
     fprintf('&%d',WideNn(i));
    for j = 1:size(K,2)
        if j == 1
            fprintf('&%d\t&%.2e\t',K(i,j),Err(i,j));
        elseif j == size(K,2)
             fprintf('&%d\t&%.2e\\\\\n',K(i,j),Err(i,j));
            
        else
           fprintf('&%d\t&%.2e\t',K(i,j),Err(i,j));
        end
    end
end