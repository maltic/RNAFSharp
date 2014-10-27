module MultiLoop

type Parameters = 
    {
        // multiloop params
        eDangleA : float;
        eDangleB : float;
        eDangleC : float;
        eDangleValues : float[]; 

        iDangleA : float;
        iDangleB : float;
        iDangleC : float;
        iDangleValues : float[];

        // multiloop stem bonuses
        eStemA : float;
        eStemB : float;
        eStemC : float;

        iStemA : float;
        iStemB : float;
        iStemC : float;
    }
    static member Random() = 
        let r = System.Random()
        let rand _ = r.NextDouble()
        let b = RNAPrimary.BASES
        {
            eDangleA = rand();
            eDangleB = rand();
            eDangleC = rand();
            eDangleValues = Array.init b rand;
            iDangleA = rand();
            iDangleB = rand();
            iDangleC = rand();
            iDangleValues = Array.init b rand;
            eStemA = rand();
            eStemB = rand();
            eStemC = rand();
            iStemA = rand();
            iStemB = rand();
            iStemC = rand();
        }

let scoreInternalDangle p (rna:RNAPrimary.Base[]) i j =
    ModelFunctions.rnaRange rna i j 
    |> ModelFunctions.jacobsonStockmayer p.iDangleA p.iDangleB p.iDangleC p.iDangleValues

let scoreExternalDangle p (rna:RNAPrimary.Base[]) i j =
    ModelFunctions.rnaRange rna i j 
    |> ModelFunctions.jacobsonStockmayer p.eDangleA p.eDangleB p.eDangleC p.eDangleValues

let scoreInternalStem p (rna:RNAPrimary.Base[]) i j =
    ModelFunctions.jacobsonStockmayerEq p.iStemA p.iStemB p.iStemC (float (j-i))

let scoreExternalStem p (rna:RNAPrimary.Base[]) i j =
    ModelFunctions.jacobsonStockmayerEq p.eStemA p.eStemB p.eStemC (float (j-i))


let mutate rate range (r:System.Random) p = 
    let maybeMut = GA.maybeMutate rate range r
    let b = RNAPrimary.BASES
    {
        eDangleA = maybeMut p.eDangleA;
        eDangleB = maybeMut p.eDangleB;
        eDangleC = maybeMut p.eDangleC;
        eDangleValues = Array.map maybeMut p.eDangleValues;
        iDangleA = maybeMut p.iDangleA;
        iDangleB = maybeMut p.iDangleB;
        iDangleC = maybeMut p.iDangleC;
        iDangleValues = Array.map maybeMut p.iDangleValues;
        eStemA = maybeMut p.eStemA;
        eStemB = maybeMut p.eStemB;
        eStemC = maybeMut p.eStemC;
        iStemA = maybeMut p.iStemA;
        iStemB = maybeMut p.iStemB;
        iStemC = maybeMut p.iStemC;
    }

let breed crossOverRate (r:System.Random) a b = 
    let maybeCross = GA.maybeTake crossOverRate r
    let bas = RNAPrimary.BASES
    {
        eDangleA = maybeCross a.eDangleA b.eDangleA;
        eDangleB = maybeCross a.eDangleB b.eDangleB;
        eDangleC = maybeCross a.eDangleC b.eDangleC;
        eDangleValues = Array.init bas (fun i -> maybeCross a.eDangleValues.[i] b.eDangleValues.[i]);
        iDangleA = maybeCross a.iDangleA b.iDangleA;
        iDangleB = maybeCross a.iDangleB b.iDangleB;
        iDangleC = maybeCross a.iDangleC b.iDangleC;
        iDangleValues = Array.init bas (fun i -> maybeCross a.iDangleValues.[i] b.iDangleValues.[i]);
        eStemA = maybeCross a.eStemA b.eStemA;
        eStemB = maybeCross a.eStemB b.eStemB;
        eStemC = maybeCross a.eStemC b.eStemC;
        iStemA = maybeCross a.iStemA b.iStemA;
        iStemB = maybeCross a.iStemB b.iStemB;
        iStemC = maybeCross a.iStemC b.iStemC;
    }