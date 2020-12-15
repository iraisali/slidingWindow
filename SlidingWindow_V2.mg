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

// ------ Conversion jacobien vers affine

function JacobianToAffine(J)
	if J eq [1,1,0] then return E!0;
	else 
		a:= J[1] /(J[3]^2);
		b:= J[2] /(J[3]^3);
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
	A:= 4*P[1]*P[2]^2;
	B:= 3*P[1]^2+a*P[3]^4;
	//Q:=[-2*A+B^2,-8*P[2]^4+B*(A-(-2*A+B^2)),2*P[2]*P[3]];
	X:= -2*A+B^2;
	Q:=[X,-8*P[2]^4+B*(A-X),2*P[2]*P[3]];
	return(Q);
end function;


// ------ Addition jacobienne
      //Entree :
                //P point
                //Q point
      //Sortie : 
                //P+Q coord jac
addJacobian:=function(P,Q) 
	P:=[P[1],P[2],P[3]];
	Q:=[Q[1],Q[2],Q[3]];
	if P eq [1,1,0] then return(Q);
	else
		if Q eq [1,1,0] then return(P);
		else
			A:= P[1]*Q[3]^2;
			B:= Q[1]*P[3]^2;
			C:= P[2]*Q[3]^3;
			D:= Q[2]*P[3]^3;
			E:= B-A;
			F:= D-C;
			X:= -E^3-2*A*E^2+F^2;
			return([X,-C*E^3+F*(A*E^2-X),P[3]*Q[3]*E]);
		end if;
	end if;
end function;


// ------ Double & multiply & addjacobien
    //En entree : //
              //n nombre de doublement
              //a coefficient de la courbe
              //P point
    //En sortie : 
              //(2n+1)*P
doubleMultiplyAdd:= function(P,n)
	y:=[P[1],P[2],P[3]];
	if n eq 0 then 
		//return([1,1,0]); //neutre
		return y;
	else
		for i:=1 to n do
			y := addJacobian(y,doubleJacobian(P));
		end for;
		return(y);
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

function PrecomputesSW(P,k)
    return [doubleMultiplyAdd(P,i) :  i in [0..2^(k-1) -1 ]];
    //return [doubleMultiplyAdd(P,i) :  i in [1..2^(k-1) -1 ]];
end function;
 

// ------------------------------------- SLIDING WINDOW -----------------------------------

function SlidingWindow(P,n,k)
	prc:=PrecomputesSW(P,k);
	L:=IntegerToSequence(n,2); //liste des coefs de n dans sa décompo en base 2
	y := [1,1,0]; //neutre
	i := #L;
	while i ge 1 do
		if L[i] eq 0 then 
			y := doubleJacobian(y); 
			i := i-1;
		else
			//s := Max(i-k+1,0);
			s := Max(i - k +1,1);
			while L[s] eq 0 do
				s := s+1;
			end while;
			for h:=1 to i-s + 1 do
				y := doubleJacobian(y);
			end for;
			u := SequenceToInteger([L[j] : j in [s..i] ] ,2);
			//y := addJacobian(y, prc[(u-1) div 2 + 1]);
			y := addJacobian(y, prc[(u+1) div 2]); //Je ne suis pas sur que ca me donne le bon
			i := s-1;
		end if;
	end while;
	return y;
end function;


function question4(n,A)
	y:= SlidingWindow(A,n,3);
	return JacobianToAffine(y);
end function;


// ----------------------------------- TESTS et TEMPS DE CALCULS --------------------------------------

B:=Random(E); n:=Random(r); n*B eq question4(n,B);
A:=G; time for i:=1 to 100 do n:=Random(r);A:=question4(n,A); end for;