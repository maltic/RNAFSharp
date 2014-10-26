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