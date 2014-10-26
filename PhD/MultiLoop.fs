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