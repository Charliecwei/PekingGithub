function q_x = QX(k,T,Dx,state)

switch state
    case 'forecast'
        q_x =  -k.*Dxscheme(T,Dx,'backward');
    case 'correct'
        q_x =  -k.*Dxscheme(T,Dx,'forward');
end

end