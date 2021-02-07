// ---------------------------------------- Paramètres ----------------------------------------

p:= 6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151;
r:= 6864797660130609714981900799081393217269435300143305409394463459185543183397655394245057746333217197532963996371363321113864768612440380340372808892707005449;
Fp:=GF(p);
b:= Fp!1093849038073734274511112390766805569936207598951683748994586394495953116150735016013708737573759623248592132296706313309438452531591012912142327488478985984;
Gx := Fp!2661740802050217063228768716723360960729859168756973147706671368418802944996427808491545080627771902352094241225065558662157113545570916814161637315895999846;
Gy := Fp!3757180025770020463545507224491183603594455134769762486694567779615544477440556316691234405012945539562144444537289428522585666729196580810124344277578376784;
a := Fp!-3;
E :=EllipticCurve([a,b]);
G:=E![Gx,Gy];

// ---------------------------------------- FONCTIONS ------------------------------------------
nu := [1,1,0];
// ------ Conversion jacobien vers affine

function JacobianToAffine(J)
	if J eq [1,1,0] then return E!0;
	else
		Z2:=J[3]^2;
		//Z3:=1/(J[3]^3);
		a:= J[1] /Z2;
		b:= J[2] /(J[3]*Z2);
		//a:=J[1]*Z3*J[3];
		//b:=J[2]*Z3;
		return E![a,b];
	end if;
end function;

// ----- Doublement Jacobien
    //En entree : //
            //P est un point en coordonees jacobiennes
            //a coefficient a de la courbe
    //En sortie : //
            //2*P avec les formule jacobiennes

doubleJacobian:= function(P)
	X:=P[1];
	Y:=P[2];
	Z:=P[3];
	YY := Y*Y;
	XYY:= X*YY;
	A:=XYY+XYY+XYY+XYY;
	XXZ4:=X*X-Z^4;
	B:=XXZ4+XXZ4+XXZ4;
	X2:=-A-A+B*B;
	Y2:=-8*YY*YY+B*(A-X2);
	YZ:=Y*Z;
    Z2:=YZ+YZ;
	return([X2,Y2,Z2]);
end function;


// ------ Addition jacobienne
      //Entree :
                //P point
                //Q point
      //Sortie : 
                //P+Q coord jac
addJacobian:=function(P,Q)
	if P eq nu then return(Q);
	elif Q eq nu then return(P);
	else
		Qz2:=Q[3]*Q[3];
		Pz2:=P[3]*P[3];
		A:= P[1]*Qz2;
		C:= -P[2]*Q[3]*Qz2;
		E:= Q[1]*Pz2-A;
		F:= Q[2]*P[3]*Pz2+C;
		e2:=E*E;
		e3:=e2*E;
		Ae2:=A*e2;
		X:=-e3-Ae2-Ae2+F*F;
		//return ([X,-C*e3+F*(A*e2-X),P[3]*Q[3]*E]);
		return ([X,C*e3+F*(Ae2-X),P[3]*Q[3]*E]);
	end if;
end function;


// --------------------------------------- PRECALCULS -------------------------------------
//Precomputes
    //En entree :
			// E une courbe elliptique
			//a coeff
			//P point
			//k taille de la fenetre
	//En sortie : //
			//liste de (2i+1)*P pour i variant de 0 a 2^(k-1)-1

precomputesSW := function(P,k)
	prc:=[P];
	S:=doubleJacobian(P);
	for i:=1 to 2^(k-1) -1 do
		prc:=Append(prc,addJacobian(prc[i],S));
	end for;
	return prc;
end function;

// ------------------------------------- SLIDING WINDOW -----------------------------------

function SlidingWindow(P,n,k)
	P:=[P[1],P[2],P[3]];
	L:=IntegerToSequence(n,2); 
	y := nu;
	i := #L;
	prc:=precomputesSW(P,k);
	while i ge 1 do
		if L[i] eq 0 then 
			y := doubleJacobian(y); 
			i := i-1;
		else
			s := Max(i - k +1,1);
			while L[s] eq 0 do
				s := s+1;
			end while;
			for h:=1 to i-s + 1 do
				y := doubleJacobian(y);
			end for;
			u := ShiftRight(ModByPowerOf2(n,i),s-1);
			y := addJacobian(y, prc[(u-1) div 2 + 1]);
			i := s-1;
		end if;
	end while;
	return y;
end function;

function question4(n,A)
	return JacobianToAffine(SlidingWindow(A,n,5));
end function;


// ----------------------------------- TESTS et TEMPS DE CALCULS --------------------------------------

B:=Random(E); n:=Random(r); n*B eq question4(n,B);
A:=G; time for i:=1 to 100 do n:=Random(r);A:=question4(n,A); end for;