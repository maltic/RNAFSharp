module Bulge

type Parameters = 
    {
        a : float;
        b : float;
        c : float;
        values : float[];
        single : float[];
    }
    static member Random() = 
        let r = System.Random()
        let rand _ = r.NextDouble()
        let b = RNAPrimary.BASES
        {
            a = rand();
            b = rand();
            c = rand();
            values = Array.init b (fun _ -> rand());
            single = Array.init b (fun _ -> rand())
        }

let score p (rna:RNAPrimary.Base[]) i j = 
    if i = j then p.single.[int rna.[i]]
    else
       ModelFunctions.jacobsonStockmayer p.a p.b p.c p.values (ModelFunctions.rnaRange rna i j)