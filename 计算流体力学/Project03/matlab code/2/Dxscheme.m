function Ds = Dxscheme(s,Dx,var)
%%X direction difference for s
Ds = zeros(size(s));
Ds(1,:) = (s(2,:)-s(1,:))/Dx;
Ds(end,:) = (s(end,:)-s(end-1,:))/Dx;
switch var
    case 'forward'
        Ds(2:end-1,:) = (s(3:end,:)-s(2:end-1,:))/Dx;
    case 'backward'
        Ds(2:end-1,:) = (s(2:end-1,:)-s(1:end-2,:))/Dx;
    case 'center'
        Ds(2:end-1,:) = (s(3:end,:)-s(1:end-2,:))/(2*Dx);
end
end