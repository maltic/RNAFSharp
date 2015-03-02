/// An energy model based on a simplified version of Tuner99.
module TurnerSimple

/// Energy parameters for stack/helix structures.
type Stack = 
    {
        init : float;
        surfaces : float[,];
    }
    /// Applies the energy model to score a stack i,j and i+1,j-1 in rna.
    member this.score (rna:RNAPrimary.Base[]) i j = 
        this.init + this.surfaces.[
            RNASecondary.basesToPair (rna.[i]) (rna.[j]) |> int,
            RNASecondary.basesToPair (rna.[i+1]) (rna.[j-1]) |> int
        ]

/// Energy parameters for hairpin loops.
type Hairpin = 
    {
        a : float;
        b : float;
        terminals : float[,,];
    }
    /// Returns the energy contribution of a hairpin loop closed by a bond i,j in a primary sequence rna.
    member this.score (rna:RNAPrimary.Base[]) i j = 
        this.terminals.[
            RNASecondary.basesToPair (rna.[i]) (rna.[j]) |> int,
             int rna.[i+1], int rna.[j-1]
        ]
        + this.a + this.b * log (j-i+1 |> float)

/// Energy parameters for bulge loops.
type Bulge = 
    {
        a : float;
        b : float;
    }
    /// Returns the energy contribution of a bulge loop of size sz.
    member this.score sz = 
        this.a + this.b * log (float sz)

/// Energy parameters for internal loops.
type Internal =
    {
        a : float;
        asymmetry : float;
        AUGUclosure : float;
        GAterminal : float;
        UUterminal : float;
    }
    /// Returns the energy contribution of an internal loop closed exteriorly by i,j and interiorly by k,j.
    member this.score (rna:RNAPrimary.Base[]) i k l j = 
        let lsz, rsz = k-i-1, j-l-1 // left and right size of internal loop
        let scoreIndexes pairs score = 
            let sorted = [| for (a,b) in pairs do yield Array.sort [|a;b|] |]
            fun a b -> 
                let s = Array.sort [|rna.[a];rna.[b]|]
                if Array.exists ((=) s) sorted then score else 0.0
        let a, u, g = RNAPrimary.Base.A, RNAPrimary.Base.U, RNAPrimary.Base.G
        ModelFunctions.jacobsonStockmayerEq 0.0 this.a 0.0 (lsz+rsz |> float)
        + if lsz <> rsz then this.asymmetry else 0.0
        + 
            let scorer = scoreIndexes [(a,u); (g,u)] this.AUGUclosure
            scorer i j + scorer k l
        + if lsz > 1 && rsz > 1 then 
                let scorerga, scoreruu = scoreIndexes [(g,a)] this.GAterminal, scoreIndexes [(u,u)] this.UUterminal
                scorerga (i+1) (j-1) + scoreruu (i+1) (j-1)
                + scorerga (k-1) (l+1) + scoreruu (k-1) (l+1) 
            else 0.0

/// Energy parameters for mutliloops.
type Multi = 
    {
        a : float;
        b : float;
        c : float;
    }
    /// Retrusn the energy contribution of a multiloop given the number of unpaired bases and branches.
    member this.score unpaired branches = 
        this.a + this.b * unpaired + this.c * branches

/// Energy parameters for the simplified Turner model.
type Parameters = 
    {
        stack : Stack;
        hp : Hairpin;
        bulge : Bulge;
        intern : Internal;
        multi : Multi;
    }
    /// Returns an instance of Parameters whose elements are extracted from stream.
    static member fromSeq (stream : float seq) = 
        let next = 
            let enumer = stream.GetEnumerator()
            fun () -> enumer.MoveNext() |> ignore; enumer.Current
        let nBases, nPairs = RNAPrimary.bases.Length, RNASecondary.pairs.Length
        {
            stack = 
                { 
                    init = next();
                    surfaces = Array2D.init nPairs nPairs (fun _ _ -> next()) 
                };
            hp = 
                {
                    a = next();
                    b = next();
                    terminals = Array3D.init nPairs nBases nBases (fun _ _ _ -> next())
                };
            bulge = 
                {
                    a = next();
                    b = next()
                };
            intern = 
                {
                    a = next();
                    asymmetry = next();
                    AUGUclosure = next();
                    GAterminal = next();
                    UUterminal = next()
                };
            multi = 
                {
                    a = next();
                    b = next();
                    c = next()
                }
        }