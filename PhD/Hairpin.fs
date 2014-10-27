module Hairpin

type Parameters =
    {
        a : float;
        b : float;
        c : float;
        values : float[];
    }
    static member Random() = 
        let r = System.Random()
        let rand _ = r.NextDouble()
        let b = RNAPrimary.BASES
        {
            a = rand();
            b = rand();
            c = rand();
            values = Array.init b rand
        }

let score p (rna:RNAPrimary.Base[]) i j = 
    ModelFunctions.rnaRange rna i j 
        |> ModelFunctions.jacobsonStockmayer p.a p.b p.c p.values


let mutate rate range (r:System.Random) p = 
    let maybeMut = GA.maybeMutate rate range r
    let b = RNAPrimary.BASES
    {
        a = maybeMut p.a;
        b = maybeMut p.b;
        c = maybeMut p.c;
        values = Array.map maybeMut p.values;
    }

let breed crossOverRate (r:System.Random) a b = 
    let maybeCross = GA.maybeTake crossOverRate r
    let bas = RNAPrimary.BASES
    {
        a = maybeCross a.a b.a
        b = maybeCross a.b b.b;
        c = maybeCross a.c b.c;
        values = Array.init bas (fun i -> maybeCross a.values.[i] b.values.[i]);
    }