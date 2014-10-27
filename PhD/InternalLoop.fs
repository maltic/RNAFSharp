module InternalLoop

type Parameters = 
    {
        a : float;
        b : float;
        c : float;
        values : float[];
        asymmetryA : float
        asymmetryB : float;
        asymmetryC : float;
    }
    static member Random() = 
        let r = System.Random()
        let rand _ = r.NextDouble()
        let b = RNAPrimary.BASES
        {
            a = rand();
            b = rand();
            c = rand();
            values = Array.init b rand;
            asymmetryA = rand();
            asymmetryB = rand();
            asymmetryC = rand();
        }

let score p (rna:RNAPrimary.Base[]) (i:int) j k l = 
    let asymmetry = System.Math.Abs((i-j) - (k-l))
    let loopScore a b = 
        ModelFunctions.jacobsonStockmayer p.a p.b p.c p.values (ModelFunctions.rnaRange rna a b)
    loopScore i j + loopScore k l +
        if asymmetry > 0 then
            ModelFunctions.jacobsonStockmayerEq p.asymmetryA p.asymmetryB p.asymmetryC (float asymmetry)
        else
            0.0


let mutate rate range (r:System.Random) p = 
    let maybeMut = GA.maybeMutate rate range r
    let b = RNAPrimary.BASES
    {
        a = maybeMut p.a;
        b = maybeMut p.b;
        c = maybeMut p.c;
        values = Array.map maybeMut p.values;
        asymmetryA = maybeMut p.asymmetryA;
        asymmetryB = maybeMut p.asymmetryB;
        asymmetryC = maybeMut p.asymmetryC;
    }

let breed crossOverRate (r:System.Random) a b = 
    let maybeCross = GA.maybeTake crossOverRate r
    let bas = RNAPrimary.BASES
    {
        a = maybeCross a.a b.a
        b = maybeCross a.b b.b;
        c = maybeCross a.c b.c;
        values = Array.init bas (fun i -> maybeCross a.values.[i] b.values.[i]);
        asymmetryA = maybeCross a.asymmetryA b.asymmetryA;
        asymmetryB = maybeCross a.asymmetryB b.asymmetryB;
        asymmetryC = maybeCross a.asymmetryC b.asymmetryC;
    }