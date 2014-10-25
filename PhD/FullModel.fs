module FullModel

type Parameters = 
    {
        // helix params
        stacking : float[,,,];
        stackA : float;
        stackB : float;
        stackC : float;
        stackEnds : float[,];

        // hairpin params
        hpA : float;
        hpB : float;
        hpC : float;
        hpValues : float[];

        // bulge params
        bulgeA : float;
        bulgeB : float;
        bulgeC : float;
        bulgeValues : float[];
        singleBulges : float[];

        // interal loop params
        ilA : float;
        ilB : float;
        ilC : float;
        ilValues : float[];
        ilAsymmetryA : float
        ilAsymmetryB : float;
        ilAsymmetryC : float;
        

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

let jacobsonStockmayerEq a b c v = 
    a*v + b*(System.Math.Log(v, System.Math.E)) + c

let jacobsonStockmayer a b c (values:float[]) data = 
    let sum = Seq.fold (fun s t -> s+values.[t]) 0.0 data
    jacobsonStockmayerEq a b c sum

let linearEq m c v = m*v + c

let stackScore p (rna:RNAPrimary.Base[]) i j k =
    jacobsonStockmayerEq p.stackA p.stackB p.stackC (float k)
        + Seq.sum 
            (seq { for l in 0..k-2 -> 
                        p.stacking.[int rna.[i+l], int rna.[j-l], int rna.[i+l+1], int rna.[j-l-1]]})
        + p.stackEnds.[int rna.[i], int rna.[j]]
        + if k > 1 then p.stackEnds.[int rna.[i+k-1], int rna.[j-k+1]] else 0.0

let hairPinScore p (rna:RNAPrimary.Base[]) i j = 
     jacobsonStockmayer p.hpA p.hpB p.hpC p.hpValues (seq {for k in i..j -> int rna.[k]})

let bulgeScore p (rna:RNAPrimary.Base[]) i j = 
    if i = j then p.singleBulges.[int rna.[i]]
    else
       jacobsonStockmayer p.bulgeA p.bulgeB p.bulgeC p.bulgeValues (seq {for k in i..j -> int rna.[k]})

let ILScore p (rna:RNAPrimary.Base[]) (i:int) j k l = 
    let asymmetry = System.Math.Abs((i-j) - (k-l))
    let loopScore a b = 
        jacobsonStockmayer p.ilA p.ilB p.ilC p.ilValues (seq {for k in a..b -> int rna.[k]})
    loopScore i j + loopScore k l +
        if asymmetry > 0 then
            jacobsonStockmayerEq p.ilAsymmetryA p.ilAsymmetryB p.ilAsymmetryC (float asymmetry)
        else
            0.0
