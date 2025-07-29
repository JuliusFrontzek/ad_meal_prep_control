import re # for regular expression operations

# Dynamikgleichung:
str1 = """c1*(xi1 - x1)*u + a11*t1*x5 + a12*t2*x6 + a13*t3*x7 + a14*t4*x8 - c2*x1 + c3*x11; 
     c1*(xi2 - x2)*u + a21*t1*x5 + a22*t2*x6 + a23*t3*x7 + a24*t4*x8 - c2*x2 + c4*x12;
     c1*(xi3 - x3)*u - a31*t1*x5 - a32*t2*x6 + a33*t3*x7 - a34*t4*x8; 
     c1*(xi4 - x4)*u - a41*t1*x5 - a42*t2*x6 - a43*t3*x7 - a44*t4*x8;
     c1*(t6*xi5 - x5)*u - t1*x5 + a55*t5*x9; 
     c1*((1-t6)*xi5 - x6)*u - t2*x6;
     c1*(xi7 - x7)*u - t3*x7 + a75*t5*x9;
     c1*(xi8 - x8)*u - t4*x8 + a85*t5*x9; 
     c1*(xi9 - x9)*u + a91*t1*x5 + a92*t2*x6 + a93*t3*x7 + a94*t4*x8 - t5*x9; 
     c1*(xi10 - x10)*u;
     c15*x11^3 + c16*x11^2*x12 + c17*x11*x12^2 + c18*x11^2 + c19*x11*x12 + c20*x11 + c5*x1; 
     c17*x12^3 + c16*x11*x12^2 + c15*x11^2*x12 + c19*x12^2 + c18*x11*x12 + c21*x12 + c5*x2"""

str2 = re.sub(r'(?<=[actx(xi)])(\d{1,2})',r'(\1)',str1)    # put in braces one or two-digit numbers after letters a, c, t, x or xi
str3 = re.sub(r'(?<=a\()(\d)',r'\1,',str2)  # put a comma in between two digits following the letter a
str4 = str3.replace('\n','')        # remove newline escape sequence
str5 = re.sub('\n\s{2,}','',str3)   # replace newline command and more than 2 whitespaces with empty string

# Messgleichung
strMgl = """(c6*x11^2 + c7*x11*x12 + c8*x12^2 + c9*x11 + c10*x12 + c11)/24;
     c12*x11; 
     c13*x12; 
     x3; 
     1 - x4/c14; 
     1 - 1E-3*x10/(c14-x4)"""
strMgl2 = re.sub(r'(?<=[actx(xi)])(\d{1,2})',r'(\1)',strMgl)    # put in braces one or two-digit numbers after letters a, c, t, x or xi
strMgl3 = strMgl2.replace('\n','')        # remove newline escape sequence

# Messgleichungen matrix-förmig:
strMglMat1 = re.sub(r'(?<=x\()(\d{1,2})',r':,\1',strMgl3) # ersetze x(10) durch x(:,10)
strMglMat2 = re.sub(r'(?=[\*\/\^])(\S)',r'.\1',strMglMat1)  # füge vor *,/,^ einen Punkt ein

# partielle Ableitung der ODEs nach x:
str_dfdx = """-c2-c1.*u,0,0,0,a11.*t1,a12.*t2,a13.*t3,a14.*t4,0,0,c3,0;
0,-c2-c1.*u,0,0,a21.*t1,a22.*t2,a23.*t3,a24.*t4,0,0,0,c4;
0,0,-c1.*u,0,-a31.*t1,-a32.*t2,a33.*t3,-a34.*t4,0,0,0,0;
0,0,0,-c1.*u,-a41.*t1,-a42.*t2,-a43.*t3,-a44.*t4,0,0,0,0;
0,0,0,0,-t1-c1.*u,0,0,0,a55.*t5,0,0,0;
0,0,0,0,0,-t2-c1.*u,0,0,0,0,0,0;
0,0,0,0,0,0,-t3-c1.*u,0,a75.*t5,0,0,0;
0,0,0,0,0,0,0,-t4-c1.*u,a85.*t5,0,0,0;
0,0,0,0,a91.*t1,a92.*t2,a93.*t3,a94.*t4,-t5-c1.*u,0,0,0;
0,0,0,0,0,0,0,0,0,-c1.*u,0,0;
c5,0,0,0,0,0,0,0,0,0,c20+2.*c18.*x11+c19.*x12+3.*c15.*x11.^2+c17.*x12.^2+2.*c16.*x11.*x12,c19.*x11+c16.*x11.^2+2.*c17.*x11.*x12;
0,c5,0,0,0,0,0,0,0,0,c18.*x12+c16.*x12.^2+2.*c15.*x11.*x12,c21+c18.*x11+2.*c19.*x12+c15.*x11.^2+3.*c17.*x12.^2+2.*c16.*x11.*x12"""

str_dfdx2 = re.sub(r'(?<=[actx(xi)])(\d{1,2})',r'(\1)',str_dfdx)    # put in braces one or two-digit numbers after letters a, c, t, x or xi
str_dfdx3 = re.sub(r'(?<=a\()(\d)',r'\1,',str_dfdx2)  # put a comma in between two digits following the letter a
str_dfdx4 = str_dfdx3.replace('\n','')        # remove newline escape sequence

# partielle Ableitung der MGL nach x:
str_dhdx = """0,0,0,0,0,0,0,0,0,0,c9./24+(c6.*x11)./12+(c7.*x12)./24,c10./24+(c7.*x11)./24+(c8.*x12)./12;
0,0,0,0,0,0,0,0,0,0,c12,0;
0,0,0,0,0,0,0,0,0,0,0,c13;
0,0,1,0,0,0,0,0,0,0,0,0;
0,0,0,-1./c14,0,0,0,0,0,0,0,0;
0,0,0,-x10./(c14-x4).^2,0,0,0,0,0,-1./(c14-x4),0,0"""

str_dhdx2 = re.sub(r'(?<=[actx(xi)])(\d{1,2})',r'(\1)',str_dhdx)    # put in braces one or two-digit numbers after letters a, c, t, x or xi
str_dhdx3 = re.sub(r'(?<=a\()(\d)',r'\1,',str_dhdx2)  # put a comma in between two digits following the letter a
str_dhdx4 = str_dhdx3.replace('\n','')
a = 1   # for debugging