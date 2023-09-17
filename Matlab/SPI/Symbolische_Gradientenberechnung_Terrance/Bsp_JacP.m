%created automatically with SPI_JacGen_Wilms
%    02-Nov-2022 16:30:04

JacP = [(mS.*mumax)./(V.*(KS+mS./V)) (mX.*mumax)./(V.*(KS+mS./V))-(mS.*mX.*mumax)./(V.^2.*(KS+mS./V).^2) (mS.^2.*mX.*mumax)./(V.^3.*(KS+mS./V).^2)-(mS.*mX.*mumax)./(V.^2.*(KS+mS./V))
-(mS.*mumax)./(V.*YXS.*(KS+mS./V)) (mS.*mX.*mumax)./(V.^2.*YXS.*(KS+mS./V).^2)-(mX.*mumax)./(V.*YXS.*(KS+mS./V)) (mS.*mX.*mumax)./(V.^2.*YXS.*(KS+mS./V))-(mS.^2.*mX.*mumax)./(V.^3.*YXS.*(KS+mS./V).^2)
0 0 0];