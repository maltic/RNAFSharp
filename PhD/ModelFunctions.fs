﻿module ModelFunctions

let jacobsonStockmayerEq a b c v = 
    a*v + b*(System.Math.Log(v, System.Math.E)) + c

let jacobsonStockmayer a b c (values:float[]) data = 
    let sum = Seq.fold (fun s t -> s+values.[t]) 0.0 data
    jacobsonStockmayerEq a b c sum

let rnaRange (rna:RNAPrimary.Base[]) i j = seq {for k in i..j -> int rna.[k]}

let linearEq m c v = m*v + c