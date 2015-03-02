/// An energy model based on a simplified version of Tuner99.
module TurnerSimple

/// Energy parameters for stack/helix structures.
type Stack = 
    {
        init : float;
        surfaces : float[,];
    }
    member this.scoreSurface (rna:RNAPrimary.Base[]) i j =
        this.surfaces.[
                RNASecondary.basesToPair (rna.[i]) (rna.[j]) |> int,
                RNASecondary.basesToPair (rna.[i+1]) (rna.[j-1]) |> int
            ]
    /// Applies the energy model to score a stack starting at i,j in rna, with k pairs.
    member this.score (rna:RNAPrimary.Base[]) i j k = 
        if k = 1 then this.init
        else
            this.scoreSurface rna i j + this.score rna (i+1) (j-1) (k-1)

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
        this.a * log (lsz+rsz |> float)
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

open RNASecondary
/// Energy parameters for the simplified Turner model.
type Parameters = 
    {
        stack : Stack;
        hp : Hairpin;
        bulge : Bulge;
        intern : Internal;
        multi : Multi;
    }

    static member toSeq p =
        let enumerate2D arr = 
            seq {
                for i in 0..Array2D.length1 arr - 1 do
                    for j in 0..Array2D.length2 arr-1 do
                        yield arr.[i,j]
            }
        let enumerate3D arr = 
            seq {
                for i in 0..Array3D.length1 arr - 1 do
                    for j in 0..Array3D.length2 arr-1 do
                        for k in 0..Array3D.length3 arr-1 do
                            yield arr.[i,j,k]
            }
        seq {
            yield p.stack.init
            yield! enumerate2D p.stack.surfaces
            yield p.hp.a
            yield p.hp.b
            yield! enumerate3D p.hp.terminals
            yield p.bulge.a
            yield p.bulge.b
            yield p.intern.a
            yield p.intern.asymmetry
            yield p.intern.AUGUclosure
            yield p.intern.GAterminal
            yield p.intern.UUterminal
            yield p.multi.a
            yield p.multi.b
            yield p.multi.c
        }

    /// Returns an instance of Parameters whose elements are extracted from stream.
    static member ofSeq (stream : float seq) = 
        let next = 
            let enumer = stream.GetEnumerator()
            fun () -> enumer.MoveNext() |> ignore; enumer.Current
        let nBases, nPairs = RNAPrimary.bases.Length, pairs.Length
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

    member p.score rna loop =
        let rec scoreInternalSurface = function
            | Reducible(children, i, j, sz) ->
                p.stack.score rna i j sz + 
                    match children with
                    | [Irreducible(_); (Reducible(_, k, l, _) as r); Irreducible(_)] -> 
                        p.intern.score rna i k l j + scoreInternalSurface r
                    | [Irreducible(i, j); (Reducible(_) as r)] -> 
                        p.bulge.score (j-i+1) + scoreInternalSurface r
                    | [(Reducible(_) as r); Irreducible(i, j)] -> 
                        p.bulge.score (j-i+1) + scoreInternalSurface r
                    | [Irreducible(_)] -> p.hp.score rna i j
                    | [] -> 0.0 // stack without hp loop
                    | l -> List.fold (fun s t -> match t with 
                                                    | Irreducible(i, j) -> p.multi.b * float (j-i+1) + s
                                                    | Reducible(_, i, j, _) as r -> 
                                                        p.multi.c + 
                                                            scoreInternalSurface r + s) 0.0 l + p.multi.a
            | Irreducible(_) -> failwith "Impossible"
        List.fold (fun s t -> match t with 
                                | Irreducible(i, j) -> p.multi.b * float (j-i+1) + s
                                | Reducible(_, i, j, _) as r -> 
                                    p.multi.c + 
                                        scoreInternalSurface r + s) 0.0 loop + p.multi.a