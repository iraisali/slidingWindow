 // Finite fields

    function PrecomputationsSW(x,k) // x est un élément d'un corps fini et k est la taille de la fenêtre.
        return [x^(2*i+1) :  i in [0..2^(k-1) -1 ]];
    end function;



    function SlidingWindow(x,n,k,prc)  // On calcule x^n avec l'algo de la fenêtre glissante, la taille de la fenêtre étant k et prc est la liste des précalculs.
        L:=IntegerToSequence(n,2);
        y := 1;
        i := #L;
        while i ge 1 do
            if L[i] eq 0 then 
                y := y^2 ; 
                i := i-1;
            else
                s := Max(i - k +1,1);
                while L[s] eq 0 do
                    s := s+1;
                end while;
                for h:=1 to i-s + 1 do
                    y := y^2;
                end for;
                u := SequenceToInteger([L[j] : j in [s..i] ] ,2);
                y := y* prc[(u-1) div 2 + 1 ];
                i := s-1;
            end if;
        end while;
        return y;
    end function;
