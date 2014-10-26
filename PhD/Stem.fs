module Stem

type Parameters =
    {
        stack : float[,,,];
        a : float;
        b : float;
        c : float;
        caps : float[,];
    }
    static member Random() = 
        let r = System.Random()
        let rand _ = r.NextDouble()
        let b = RNAPrimary.BASES
        {
            stack = Array4D.init b b b b (fun _ _ _ _ -> rand ())
            a = rand();
            b = rand();
            c = rand();
            caps = Array2D.init b b (fun _ _ -> rand())
        }


let score p (rna:RNAPrimary.Base[]) i j k =
    ModelFunctions.jacobsonStockmayerEq p.a p.b p.c (float k)
        + Seq.sum 
            (seq { for l in 0..k-2 -> 
                        p.stack.[int rna.[i+l], int rna.[j-l], int rna.[i+l+1], int rna.[j-l-1]]})
        + p.caps.[int rna.[i], int rna.[j]]
        + if k > 1 then p.caps.[int rna.[i+k-1], int rna.[j-k+1]] else 0.0

let mutate rate range (r:System.Random) p = 
    let maybeMut = GA.maybeMutate rate range r
    let b = RNAPrimary.BASES
    {
        stack = Array4D.init b b b b (fun a b c d -> maybeMut p.stack.[a,b,c,d]);
        a = maybeMut p.a;
        b = maybeMut p.b;
        c = maybeMut p.c;
        caps = Array2D.map maybeMut p.caps
    }

let breed crossOverRate (r:System.Random) a b = 
    let maybeCross = GA.maybeTake crossOverRate r
    let bas = RNAPrimary.BASES
    {
        stack = Array4D.init bas bas bas bas (fun i j k l -> maybeCross a.stack.[i,j,k,l] b.stack.[i,j,k,l]);
        a = maybeCross a.a b.a
        b = maybeCross a.b b.b;
        c = maybeCross a.c b.c;
        caps = Array2D.init bas bas (fun i j -> maybeCross a.caps.[i,j] b.caps.[i,j])
    }