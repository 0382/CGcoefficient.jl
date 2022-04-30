var documenterSearchIndex = {"docs":
[{"location":"formula/#Formula","page":"Formula","title":"Formula","text":"","category":"section"},{"location":"formula/#CG-coefficient","page":"Formula","title":"CG coefficient","text":"","category":"section"},{"location":"formula/","page":"Formula","title":"Formula","text":"Ref [1], P240, Section 8.2.4, Formula (20).","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"C_j_1m_1j_2m_2^j_3m_3 = delta_m_3 m_1+m_2\nleft\n    dfracbeginpmatrix2j_1  J-2j_2endpmatrixbeginpmatrix2j_2J-2j_3endpmatrixbeginpmatrixJ+1J-2j_3endpmatrixbeginpmatrix2j_1j_1-m_1endpmatrixbeginpmatrix2j_2  j_2 - m_2endpmatrixbeginpmatrix2j_3  j_3-m_3endpmatrix\nright^12 \ntimes sum_z (-1)^zbeginpmatrixJ-2j_3  zendpmatrixbeginpmatrixJ-2j_2  j_1-m_1-zendpmatrixbeginpmatrixJ-2j_1  j_2 + m_2 - zendpmatrix","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"where, J = j_1+j_2+j_3. It is already combination of binomials.","category":"page"},{"location":"formula/#j-symbol","page":"Formula","title":"3j symbol","text":"","category":"section"},{"location":"formula/","page":"Formula","title":"Formula","text":"Ref [1], P236, Section 8.1.2, Formula (11).","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"beginpmatrixj_1  j_2  j_3  m_1  m_2  m_3endpmatrix = (-1)^j_3+m_3+2j_1dfrac1sqrt2j_3+1C_j_1-m_1j_2-m_2^j_3m_3","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"This package use the CG coefficient above to calculate 3j symbol.","category":"page"},{"location":"formula/#j-symbol-2","page":"Formula","title":"6j symbol","text":"","category":"section"},{"location":"formula/","page":"Formula","title":"Formula","text":"Ref [1], P293, Section 9.2.1, Formula (1).","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"beginBmatrixj_1  j_2  j_3  j_4  j_5  j_6endBmatrix = Delta(j_1j_2j_3)Delta(j_4j_5j_3)Delta(j_1j_5j_6)Delta(j_4j_2j_6) \ntimes sumlimits_xdfrac(-1)^x(x+1)(x-j_123)(x-j_453)(x-j_156)(x-j_426)(j_1245-x)(j_1346-x)(j_2356-x)","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"Here, I use j_123 equiv j_1+j_2+j_3 for simplicity. The symbol Delta(abc) is defined as","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"Delta(abc) = leftdfrac(a+b-c)(a-b+c)(-a+b+c)(a+b+c+1)right^frac12","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"We can find that","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"beginpmatrixj_1+j_2-j_3  x - j_453endpmatrix = dfrac(j_1+j_2-j_3)(x-j_453)(j_1245-x) \nbeginpmatrixj_1-j_2+j_3  x - j_426endpmatrix = dfrac(j_1-j_2+j_3)(x-j_426)(j_1346-x) \nbeginpmatrixj_2+j_3-j_1  x - j_156endpmatrix = dfrac(j_2+j_3-j_1)(x-j_156)(j_2356-x) \nbeginpmatrixx+1  j_123+1endpmatrix = dfrac(x+1)(x-j_123)(j_123+1)","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"So, we have","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"beginBmatrixj_1  j_2  j_3  j_4  j_5  j_6endBmatrix = dfracDelta(j_4j_5j_3)Delta(j_1j_5j_6)Delta(j_4j_2j_6)Delta(j_1j_2j_3) \ntimes sumlimits_x (-1)^x beginpmatrixx+1  j_123+1endpmatrix beginpmatrixj_1+j_2-j_3  x - j_453endpmatrix beginpmatrixj_1-j_2+j_3  x - j_426endpmatrix beginpmatrixj_2+j_3-j_1  x - j_156endpmatrix","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"Rewrite Delta(abc) with binomials,","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"Delta(abc) = leftdfrac1beginpmatrixa+b+c+1  2a + 1endpmatrix beginpmatrix2a  a + b - cendpmatrix(2a+1)right^frac12","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"So","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"dfracDelta(j_4j_5j_3)Delta(j_1j_5j_6)Delta(j_4j_2j_6)Delta(j_1j_2j_3) \n= dfrac12j_4+1 leftdfracbeginpmatrixj_123+1  2j_1+1endpmatrixbeginpmatrix2j_1  j_1 + j_2 - j_3endpmatrixbeginpmatrixj_156+1  2j_1+1endpmatrixbeginpmatrix2j_1  j_1+j_5-j_6endpmatrixbeginpmatrixj_453+12j_4+1endpmatrixbeginpmatrix2j_4j_4+j_5-j_3endpmatrixbeginpmatrixj_426+1  2j_4+1endpmatrixbeginpmatrix2j_4  j_4+j_2-j_6endpmatrixright^frac12","category":"page"},{"location":"formula/#Racha-coefficient","page":"Formula","title":"Racha coefficient","text":"","category":"section"},{"location":"formula/","page":"Formula","title":"Formula","text":"Ref [1], P291, Section 9.1.2, Formula (11)","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"W(j_1j_2j_3j_4 j_5j_6) = (-1)^j_1+j_2+j_3+j_4 beginBmatrix\nj_1  j_2  j_5 \nj_4  j_3  j_6\nendBmatrix","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"This package uses the 6j symbol above to calculate Racha coefficient. ","category":"page"},{"location":"formula/#j-symbol-3","page":"Formula","title":"9j symbol","text":"","category":"section"},{"location":"formula/","page":"Formula","title":"Formula","text":"Ref [1], P340, Section 10.2.4, Formula (20)","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"beginBmatrixj_1  j_2  j_3  j_4  j_5  j_6  j_7  j_8  j_9endBmatrix = sumlimits_t(-1)^2t(2t+1)beginBmatrixj_1  j_2  j_3  j_6  j_9  tendBmatrix beginBmatrixj_4  j_5  j_6  j_2  t  j_8endBmatrix beginBmatrixj_7  j_8  j_9  t  j_1  j_4endBmatrix","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"Use the 6j symbol result above, we get","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"dfracDelta(j_1j_9t)Delta(j_6j_9j_3)Delta(j_6j_2t)Delta(j_1j_2j_3) dfracDelta(j_4tj_8)Delta(j_2tj_6)Delta(j_2j_5j_8)Delta(j_4j_5j_6) dfracDelta(j_7j_1j_4)Delta(tj_1j_9)Delta(tj_8j_4)Delta(j_7j_8j_9) \n = dfracDelta(j_3j_6j_9)Delta(j_2j_5j_8)Delta(j_1j_4j_7)Delta(j_1j_2j_3)Delta(j_4j_5j_6)Delta(j_7j_8j_9) Delta^2(j_1j_9t)Delta^2(j_2j_6t)Delta^2(j_4j_8t)","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"Define","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"P_0 equiv dfracDelta(j_3j_6j_9)Delta(j_2j_5j_8)Delta(j_1j_4j_7)Delta(j_1j_2j_3)Delta(j_4j_5j_6)Delta(j_7j_8j_9) \n= leftdfracbeginpmatrixj_123 + 1  2j_1+1endpmatrixbeginpmatrix2j_1j_1+j_2-j_3endpmatrixbeginpmatrixj_456+12j_5+1endpmatrixbeginpmatrix2j_5  j_4+j_5-j_6endpmatrixbeginpmatrixj_789+1  2j_9+1endpmatrixbeginpmatrix2j_9  j_7 + j_9 - j_8endpmatrixbeginpmatrixj_147 + 1 2j_1 + 1endpmatrixbeginpmatrix2j_1  j_1+j_4-j_7endpmatrixbeginpmatrixj_258+1  2j_5+1endpmatrixbeginpmatrix2j_5  j_2+j_5 - j_8endpmatrixbeginpmatrixj_369+1  2j_9+1endpmatrixbeginpmatrix2j_9  j_3+j_9 - j_6endpmatrixright^12","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"and then define","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"P(t) equiv (-1)^2t(2t+1)Delta^2(j_1j_9t)Delta^2(j_2j_6t)Delta^2(j_4j_8t)\n = dfrac(-1)^2t(2t+1)^2 times \n\ndfrac1beginpmatrixj_1+j_9+t+1  2t+1endpmatrixbeginpmatrix2t  j_1+t-j_9endpmatrixbeginpmatrixj_2+j_6+t+1  2t+1endpmatrixbeginpmatrix2t  j_2+t-j_6endpmatrixbeginpmatrixj_4+j_8+t  2t+1endpmatrixbeginpmatrix2t  j_4+t-j_8endpmatrix","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"A(tx) equiv (-1)^xbeginpmatrixx+1  j_123 + 1endpmatrixbeginpmatrixj_1+j_2-j_3  x - (j_6+j_9+j_3)endpmatrix beginpmatrixj_1+j_3-j_2  x - (j_6+j_2+t)endpmatrixbeginpmatrixj_2+j_3-j_1  x - (j_1 + j_9 + t)endpmatrix","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"B(ty) equiv (-1)^ybeginpmatrixy+1  j_456 + 1endpmatrixbeginpmatrixj_4+j_5-j_6  y - (j_2+t+j_6)endpmatrixbeginpmatrixj_4+j_6-j_5  y - (j_2+j_5+j_8)endpmatrixbeginpmatrixj_5+j_6-j_4  y - (j_4+t+j_8)endpmatrix","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"C(tz) equiv (-1)^z beginpmatrixz+1  j_789 + 1endpmatrixbeginpmatrixj_7+j_8-j_9  z - (t+j_1+j_9)endpmatrixbeginpmatrixj_7+j_9-j_8  z - (t+j_8+j_4)endpmatrixbeginpmatrixj_8+j_9- j_7  z - (j_7 + j_1 + j_4)endpmatrix","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"At last, we get","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"beginBmatrixj_1  j_2  j_3  j_4  j_5  j_6  j_7  j_8  j_9endBmatrix = P_0sum_t P(t) left(sum_x A(tx)right) left(sum_y B(ty)right) left(sum_z C(tz)right)","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"It deserves to be mentioned that, although the formula has 4 sums, the sum of xyz are decoupled. So we can do the three for loops respectively, which means the depth of for loop is not 4 but 2.","category":"page"},{"location":"formula/#Estimate-the-capacity","page":"Formula","title":"Estimate the capacity","text":"","category":"section"},{"location":"formula/","page":"Formula","title":"Formula","text":"Assume we are doing a calculation for a system, we will not calculate the Winger Symbols with very large angular momentum, because usually we will take some trunction. If in such truncated system, the max angular momentum is J_max, now we can estimate how many binomial coefficients we need to store to compute those Wigner Symbols.","category":"page"},{"location":"formula/#With-maximum-J_{max}","page":"Formula","title":"With maximum J_max","text":"","category":"section"},{"location":"formula/#CG-and-3j","page":"Formula","title":"CG & 3j","text":"","category":"section"},{"location":"formula/","page":"Formula","title":"Formula","text":"According to the formula, the max possible binomial is beginpmatrix J+1  J-2j_3endpmatrix, where J = j_1+j_2+j_3. So we need at least store binomials to beginpmatrix n_min  kendpmatrix where","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"n_min geq 3J_max+1","category":"page"},{"location":"formula/#j-and-Racha","page":"Formula","title":"6j & Racha","text":"","category":"section"},{"location":"formula/","page":"Formula","title":"Formula","text":"For binomials in sqrt，the maximum possible n = 3J_max+1, and for those in the summation, we have to calculate the boundary of x + 1. According to the formula","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"x leq minj_1+j_2+j_4+j_5 j_1+j_3+j_4+j_6 j_2+j_3+j_5+j_6 leq 4J_max","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"So we have to at least store to","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"n_min geq 4J_max + 1","category":"page"},{"location":"formula/#j","page":"Formula","title":"9j","text":"","category":"section"},{"location":"formula/","page":"Formula","title":"Formula","text":"Similarly, according to P_0, we need 3J_max+1. According to P(t),","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"t leq minj_2 + j_6 j_4 + j_8 j_1 + j_9 leq 2J_max","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"so we need j_1 + j_9 + t + 1 leq 4J_max + 1. According to A(t x) B(t y) C(t z),","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"x leq minj_1+j_2+j_6+j_9 j_1+j_3+j_6+t j_2+j_3+j_9+t leq 5J_max","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"so we need to at least store to","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"n_min geq 5J_max + 1","category":"page"},{"location":"formula/#With-maximum-single-particle-j_{max}","page":"Formula","title":"With maximum single particle j_max","text":"","category":"section"},{"location":"formula/","page":"Formula","title":"Formula","text":"Considering a quantum many body calculation, we can get a more optimistic estimation of the calculation capacity. In a quantum many body calculation, we often truncate single particle orbits, and in this condition we assume the maximum angular momentum of the single particle orbits is j_max.","category":"page"},{"location":"formula/#CG-and-3j-2","page":"Formula","title":"CG & 3j","text":"","category":"section"},{"location":"formula/","page":"Formula","title":"Formula","text":"In this condition, the j_1 j_2 j_3 contains two single particle angular momentum, and the rest one is a coupled two body angular momentum. So","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"n_min geq 4j_max + 1","category":"page"},{"location":"formula/#j-and-Racha-2","page":"Formula","title":"6j & Racha","text":"","category":"section"},{"location":"formula/","page":"Formula","title":"Formula","text":"Similarly","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"x leq minj_1+j_2+j_4+j_5 j_1+j_3+j_4+j_6 j_2+j_3+j_5+j_6","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"These summations of four angular momentum can be split into two kinds, all single particle angular momentum, and two single particle with two two body angular momentum. Consider the worst condition","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"n_min geq 6j_max + 1","category":"page"},{"location":"formula/#j-2","page":"Formula","title":"9j","text":"","category":"section"},{"location":"formula/","page":"Formula","title":"Formula","text":"We assume the 9j symbol contains at least one single particle angular momentum, then according to the coupling rules, the 9j symbol can be split into kinds","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"beginBmatrix\n1  1  2 \n1  1  2 \n2  2  2\nendBmatrix\nquad\nbeginBmatrix\n1  1  2 \n1  2  1 \n2  1  1\nendBmatrix","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"where 1 represents single particle and 2 for two body angular momentum. In both cases, we can deduce that","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"n_min geq 8j_max + 1","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"Reference","category":"page"},{"location":"formula/","page":"Formula","title":"Formula","text":"[1]: A. N. Moskalev D. A. Varshalovich and V. K. Khersonskii, Quantum theory of angular momentum.","category":"page"},{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/#Types","page":"API","title":"Types","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"HalfInt\nSqrtRational","category":"page"},{"location":"api/#CGcoefficient.HalfInt","page":"API","title":"CGcoefficient.HalfInt","text":"HalfInt = Union{Integer, Rational}\n\nAngular momentum quantum numbers may be half integers and integers. With HalfInt type, you can use both 2 and 3//2 as parameters. But the parameter like 5//3, who's denominator is not 2 while gives out error.\n\n\n\n\n\n","category":"type"},{"location":"api/#CGcoefficient.SqrtRational","page":"API","title":"CGcoefficient.SqrtRational","text":"SqrtRational <: Real\n\nStore exact value of a Rational's square root. It is stored in s√r format, and we do not simplify it during arithmetics. You can use the simplify function to simplify it.\n\n\n\n\n\n","category":"type"},{"location":"api/#Core-functions","page":"API","title":"Core functions","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"CG\nthreeJ\nsixJ\nnineJ\nRacah","category":"page"},{"location":"api/#CGcoefficient.CG","page":"API","title":"CGcoefficient.CG","text":"CG(j1::HalfInt, j2::HalfInt, j3::HalfInt, m1::HalfInt, m2::HalfInt, m3::HalfInt)\n\nCG coefficient langle j_1m_1 j_2m_2  j_3m_3 rangle\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.threeJ","page":"API","title":"CGcoefficient.threeJ","text":"threeJ(j1::HalfInt, j2::HalfInt, j3::HalfInt, m1::HalfInt, m2::HalfInt, m3::HalfInt)\n\nWigner 3j-symbol\n\nbeginpmatrix\nj_1  j_2  j_3 \nm_1  m_2  m_3\nendpmatrix\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.sixJ","page":"API","title":"CGcoefficient.sixJ","text":"sixJ(j1::HalfInt, j2::HalfInt, j3::HalfInt, j4::HalfInt, j5::HalfInt, j6::HalfInt)\n\nWigner 6j-symbol\n\nbeginBmatrix\nj_1  j_2  j_3 \nj_4  j_5  j_6\nendBmatrix\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.nineJ","page":"API","title":"CGcoefficient.nineJ","text":"nineJ(j1::HalfInt, j2::HalfInt, j3::HalfInt,\n      j4::HalfInt, j5::HalfInt, j6::HalfInt,\n      j7::HalfInt, j8::HalfInt, j9::HalfInt)\n\nWigner 9j-symbol\n\nbeginBmatrix\nj_1  j_2  j_3 \nj_4  j_5  j_6 \nj_7  j_8  j_9\nendBmatrix\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.Racah","page":"API","title":"CGcoefficient.Racah","text":"Racah(j1::HalfInt, j2::HalfInt, j3::HalfInt, j4::HalfInt, j5::HalfInt, j6::HalfInt)\n\nRacah coefficient\n\nW(j_1j_2j_3j_4 j_5j_6) = (-1)^j_1+j_2+j_3+j_4 beginBmatrix\nj_1  j_2  j_5 \nj_4  j_3  j_6\nendBmatrix\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"API","title":"API","text":"People often use double of angular momentum quantum number as parameters, so we can use integer as parameters. This package also offers such functions, where the d letter means double.","category":"page"},{"location":"api/","page":"API","title":"API","text":"dCG\nd3j\nd6j\nd9j\ndRacah","category":"page"},{"location":"api/#CGcoefficient.dCG","page":"API","title":"CGcoefficient.dCG","text":"dCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)\n\nCG coefficient function with double angular monentum number parameters, so that the parameters can be integer. You can calculate dCG(1, 1, 2, 1, 1, 2) to calculate the real CG(1//2, 1//2, 1, 1/2, 1//2, 1)\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.d3j","page":"API","title":"CGcoefficient.d3j","text":"d3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)\n\n3j-symbol function with double angular monentum number parameters, so that the parameters can be integer.\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.d6j","page":"API","title":"CGcoefficient.d6j","text":"d6j(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)\n\n6j-symbol function with double angular monentum number parameters, so that the parameters can be integer.\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.d9j","page":"API","title":"CGcoefficient.d9j","text":"d9j(dj1::Integer, dj2::Integer, dj3::Integer,\n    dj4::Integer, dj5::Integer, dj6::Integer,\n    dj7::Integer, dj8::Integer, dj9::Integer)\n\n9j-symbol function with double angular monentum number parameters, so that the parameters can be integer.\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.dRacah","page":"API","title":"CGcoefficient.dRacah","text":"dRacah(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)\n\nRacah coefficient function with double angular momentum parameters, so that the parameters can be integer.\n\n\n\n\n\n","category":"function"},{"location":"api/#float-version-functions","page":"API","title":"float version functions","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Float version functions is always used for numeric calculation, so the parameters of all these functions (except reserve_fbinomial) are double of the exact angular momentum quantum number.","category":"page"},{"location":"api/","page":"API","title":"API","text":"fbinomial\nfCG\nf3j\nf6j\nf9j\nfRacah\nreserve_fbinomial","category":"page"},{"location":"api/#CGcoefficient.fbinomial","page":"API","title":"CGcoefficient.fbinomial","text":"fbinomial(n::Integer, k::Integer)\n\nbinomial with Float64 return value.\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.fCG","page":"API","title":"CGcoefficient.fCG","text":"fCG(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)\n\nfloat64 and fast CG coefficient.\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.f3j","page":"API","title":"CGcoefficient.f3j","text":"f3j(dj1::Integer, dj2::Integer, dj3::Integer, dm1::Integer, dm2::Integer, dm3::Integer)\n\nfloat64 and fast Wigner 3j symbol.\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.f6j","page":"API","title":"CGcoefficient.f6j","text":"f6j(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)\n\nfloat64 and fast Wigner 6j symbol.\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.f9j","page":"API","title":"CGcoefficient.f9j","text":"f9j(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer, dj7::Integer, dj8::Integer, dj9::Integer)\n\nfloat64 and fast Wigner 9j symbol.\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.fRacah","page":"API","title":"CGcoefficient.fRacah","text":"Racah(dj1::Integer, dj2::Integer, dj3::Integer, dj4::Integer, dj5::Integer, dj6::Integer)\n\nfloat64 and fast Racah coefficient.\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.reserve_fbinomial","page":"API","title":"CGcoefficient.reserve_fbinomial","text":"reserve_fbinomial(n::Integer, mode::AbstractString, rank::Integer)\n\nThis function reserves memory for fbinomial(n, k). In this code, the fbinomial function is only valid in the stored range. If you call a fbinomial function out of the range, it just gives you 0. The __init__() function stores to nmax = 67.\n\nThe parameters means\n\n Calculate range CG & 3j 6j & Racah 9j\nmeaning of type type\\rank 3 6 9\nmax angular momentum \"Jmax\" 3*Jmax+1 4*Jmax+1 5*Jmax+1\nmax two-body coupled angular momentum \"2bjmax\" 2*jmax+1 3*jmax+1 4*jmax+1\nmax binomial \"nmax\" nmax namx nmax\n\nThe \"2bjmax\" mode means your calculation only consider two-body coupling, and no three-body coupling. This mode assumes that in all these coefficients, at least one of the angular momentun is just a single particle angular momentum. With this assumption, \"2bjmax\" mode will use less memory than \"Jmax\" mode.\n\n\"Jmax\" means the global maximum angular momentum, for every parameters. It is always safe with out any assumption.\n\nThe \"nmax\" mode directly set nmax, and the rank parameter is ignored. \n\nrank = 6 means you only need to calculate CG and/or 6j symbols, you don't need to calculate 9j symbol.\n\nFor example\n\nreserve_fbinomial(21, \"Jmax\", 6)\n\nmeans you calculate CG and 6j symbols, and donot calculate 9j symbol. The maximum angular momentum in your system is 21.\n\nYou do not need to rememmber those values in the table. You just need to find the maximum angular momentum in you canculation, then call the function.\n\nThe reserve_fbinomial function is not thread safe, so you should call it before you start your calculation.\n\n\n\n\n\n","category":"function"},{"location":"api/#Some-useful-function","page":"API","title":"Some useful function","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"iphase\nis_same_parity\ncheck_jm\ncheck_couple\nexact_sqrt\nsimplify(::Integer)\nsimplify(::SqrtRational)","category":"page"},{"location":"api/#CGcoefficient.iphase","page":"API","title":"CGcoefficient.iphase","text":"iphase(n::Integer)\n\n(-1)^n\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.is_same_parity","page":"API","title":"CGcoefficient.is_same_parity","text":"is_same_parity(x::T, y::T) where {T <: Integer}\n\njudge if two integers are same odd or same even\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.check_jm","page":"API","title":"CGcoefficient.check_jm","text":"check_jm(dj::T, dm::T) where {T <: Integer}\n\ncheck if the m-quantum number if one of the components of the j-quantum number, in other words, m and j has the same parity, and abs(m) < j\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.check_couple","page":"API","title":"CGcoefficient.check_couple","text":"check_couple(dj1::T, dj2::T, dj3::T) where {T <: Integer}\n\ncheck if three angular monentum number j1, j2, j3 can couple\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.exact_sqrt","page":"API","title":"CGcoefficient.exact_sqrt","text":"exact_sqrt(r::Union{Integer, Rational})\n\nGet exact √r using SqrtRational type.\n\n\n\n\n\n","category":"function"},{"location":"api/#CGcoefficient.simplify-Tuple{Integer}","page":"API","title":"CGcoefficient.simplify","text":"simplify(n::Integer)\n\nSimplify a integer n = x * t^2 to (x, t)\n\n\n\n\n\n","category":"method"},{"location":"api/#CGcoefficient.simplify-Tuple{SqrtRational}","page":"API","title":"CGcoefficient.simplify","text":"simplify(x::SqrtRational)\n\nSimplify a SqrtRational.\n\n\n\n\n\n","category":"method"},{"location":"#Home","page":"Home","title":"Home","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A package to calculate CG-coefficient, Racha coefficient, and Wigner 3j, 6j, 9j symbols. It store the exact result with SqrtRational type.","category":"page"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package is inspired by Ref [1]. See CENS-MBPT for details.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The idea is to simplify 3nj Symbols to sum combinations of binomial coefficients. We can calculate binomial coefficients by Pascal's Triangle, and store them first. Then we calculate 3nj Symbols using the stored binomial coefficients. However, in current version, I just use the builtin binomial function, and calculation with large integer will overflow.","category":"page"},{"location":"#Install","page":"Home","title":"Install","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Just install with Julia REPL and enjoy it.","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add CGcoefficient","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"push!(LOAD_PATH, \"../../src/\") # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"using CGcoefficient\nsixJ(1,2,3,4,5,6)","category":"page"},{"location":"","page":"Home","title":"Home","text":"In a markdown enviroment, such as jupyter notebook, it will give you a latex output. For the seek of speed, we do not simplify the result during the arithmetics. You can simplify the result explicitly.","category":"page"},{"location":"","page":"Home","title":"Home","text":"simplify(sixJ(1,2,3,4,5,6))","category":"page"},{"location":"","page":"Home","title":"Home","text":"You can also do some arithmetics with the result, thus do arithmetics using the SqrtRational type.","category":"page"},{"location":"","page":"Home","title":"Home","text":"x = sixJ(1,2,3,4,5,6) * exact_sqrt(1//7) * exact_sqrt(1//13) * iphase(2+3+5+6)\nsimplify(x)","category":"page"},{"location":"","page":"Home","title":"Home","text":"In a console enviroment it will give out a text output.","category":"page"},{"location":"","page":"Home","title":"Home","text":"simplify(nineJ(1,2,3,5,4,3,6,6,0))","category":"page"},{"location":"","page":"Home","title":"Home","text":"You can also use print function to force print a text output.","category":"page"},{"location":"","page":"Home","title":"Home","text":"a = Racah(1,2,3,2,1,2)\nprint(simplify(a))","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Formula\nAPI","category":"page"},{"location":"","page":"Home","title":"Home","text":"[1]: T. Engeland and M. Hjorth-Jensen, the Oslo-FCI code. https://github.com/ManyBodyPhysics/CENS.","category":"page"}]
}
