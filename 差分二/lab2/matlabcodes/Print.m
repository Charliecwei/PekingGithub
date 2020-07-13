function Print(N,P,Errinf,Errinf_Rate,Wp2Err,Wp2Err_Rate,InteK,Hises)
  fprintf('h\t& Errinf\t &Errinf_Rate\t &ErrW12\t& ErrW12_Rate\t &ErrW22\t&ErrW22_Rate\\\\\n');
for i = 1:length(N)
     fprintf('%.2f&%d',2/N(i),InteK(i));
    for j = 1:(length(P)+1)
        if j == 1
            fprintf('&%.2e\t&%.2f\t',Errinf(i),Errinf_Rate(i));
        elseif j == (length(P)+1)
             fprintf('&%.2e\t&%.2f\\\\\n',Wp2Err(i,j-1),Wp2Err_Rate(i,j-1));
            
        else
           fprintf('&%.2e\t&%.2f\t',Wp2Err(i,j-1),Wp2Err_Rate(i,j-1));
        end
    end
end

for k = 1:length(N)
    K = (log10(Hises{k}(end))-log10(Hises{k}(1)))/length(Hises{k});
    h = 2/N(k);
    if k~=5
        subplot(3,2,k)
    else
        subplot(3,2,[5,6])
    end 
    plot(log10(Hises{k}),'-*');
    ylabel('log_{10}(err)');
    xlabel('µü´ú´ÎÊý');
    title(['h=',num2str(h),' k = ',num2str(K)])
end