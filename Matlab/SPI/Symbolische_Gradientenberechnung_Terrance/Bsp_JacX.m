%created automatically with SPI_JacGen_Wilms
%    02-Nov-2022 16:30:04

JacX = [(mS.*mX)./(V.*(KS+mS./V)),-(mS.*mX.*mumax)./(V.*(KS+mS./V).^2),0;
-(mS.*mX)./(V.*YXS.*(KS+mS./V)),(mS.*mX.*mumax)./(V.*YXS.*(KS+mS./V).^2),(mS.*mX.*mumax)./(V.*YXS.^2.*(KS+mS./V));
0,0,0];