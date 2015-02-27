module ModelFunctions

let jacobsonStockmayerEq a b c v = 
    a*v + b*(if v <= 0.0 then 0.0 else System.Math.Log(v, System.Math.E)) + c

let jacobsonStockmayer a b c (values:float[]) data = 
    jacobsonStockmayerEq a b c (Array.sum values)

let rnaRange (rna:RNAPrimary.Base[]) i j = seq {for k in i..j -> int rna.[k]}

let linearEq m c v = m*v + c