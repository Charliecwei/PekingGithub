function Dphip = getDphip3DP3N(Dlambda,lambda)

%Dla1,Dla2,Dla3,Dla4,la1,la2,la3,la4;
Dla1 = Dlambda(:,:,1);
Dla2 = Dlambda(:,:,2);
Dla3 = Dlambda(:,:,3);
Dla4 = Dlambda(:,:,4);
la1 = lambda(1);
la2 = lambda(2);
la3 = lambda(3);
la4 = lambda(4);





Dphip(:,:,1)=dingdian(la1,Dla1);
Dphip(:,:,2)=dingdian(la2,Dla2);
Dphip(:,:,3)=dingdian(la3,Dla3);
Dphip(:,:,4)=dingdian(la4,Dla4);

Dphip(:,:,5)=mian(la2,la3,la4,Dla2,Dla3,Dla4);
Dphip(:,:,6)=mian(la3,la4,la1,Dla3,Dla4,Dla1);
Dphip(:,:,7)=mian(la4,la1,la2,Dla4,Dla1,Dla2);
Dphip(:,:,8)=mian(la1,la2,la3,Dla1,Dla2,Dla3);

Dphip(:,:,[9,10])=bian(la1,la2,Dla1,Dla2);

Dphip(:,:,[11,12])=bian(la1,la3,Dla1,Dla3);

Dphip(:,:,[13,14])=bian(la1,la4,Dla1,Dla4);

Dphip(:,:,[15,16])=bian(la2,la3,Dla2,Dla3);

Dphip(:,:,[17,18])=bian(la2,la4,Dla2,Dla4);

Dphip(:,:,[19,20])=bian(la3,la4,Dla3,Dla4);

Dphip(:,:,21) = 18*(la1.*Dla1+la2.*Dla2+la3.*Dla3+la4.*Dla4)...
                           - 33*(la1.^2.*Dla1+la2.^2.*Dla2+la3.^2.*Dla3+la4.^2.*Dla4)...
                           - 72*(Dla1.*la2.*la3+la1.*Dla2.*la3+la1.*la2.*Dla3...
                                    +Dla1.*la2.*la4+la1.*Dla2.*la4+la1.*la2.*Dla4...
                                    +Dla1.*la3.*la4+la1.*Dla3.*la4+la1.*la3.*Dla4...
                                    +Dla2.*la3.*la4+la2.*Dla3.*la4+la2.*la3.*Dla4);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%以la=lambda(:,i)，Dla=Dlambda(:,:,i)算顶点i的基函数导数
function  Dphips=dingdian(la,Dla)
Dphips=0.5*Dla.*(27*power(la,2)-18*la+2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%算'1,2,3'点对应的面上的基函数导数
function Dphips=mian(la1,la2,la3,Dla1,Dla2,Dla3)
   Dphips=27*(Dla1.*la2.*la3+Dla2.*la3.*la1+Dla3.*la1.*la2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%算'1,2'边上的基函数导数
    function Dphips=bian(la1,la2,Dla1,Dla2)
        Dphips(:,:,1)=0.5*Dla2.*(27*power(la1,2)-9*la1)+...
                         0.5*la2.*(54*la1-9).*Dla1;
        Dphips(:,:,2)=0.5*Dla1.*(27*power(la2,2)-9*la2)+...
                         0.5*la1.*(54*la2-9).*Dla2;
    end
end
