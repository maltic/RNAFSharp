module MFE

open RNASecondary
open System

let zuker (rna:RNAPrimary.Base[]) model = 
    let dpW = Array2D.create rna.Length rna.Length Double.MinValue
    let dpV = Array2D.create rna.Length rna.Length Double.MinValue
    let dpIM = Array2D.create rna.Length rna.Length Double.MinValue
    let dpEM = Array.create rna.Length Double.MinValue
    let IL_MAX_SIZE = 30
    let rec W i j = // i and j begin a loop (i-1 and j+1 are bonded)
        if i > j then Double.MinValue
        elif dpW.[i,j] > Double.MinValue then dpW.[i,j]
        else
            let mutable max = model.hairpin i j // hairpin
            for k in i..(j-1) do // k is the start of the succeeding stem, bulge from i, k-1
                let tmp = model.lbulge i (k-1) + V k j
                if tmp > max then max <- tmp
            for k in j..(-1)..i+1 do // k is the end of the precding stem, bulge from k+1, j
                let tmp = model.rbulge (k+1) j + V i k
                if tmp > max then max <- tmp
            // internal loops
            for k in i..Math.Min(j-3, i+IL_MAX_SIZE) do // left loop i, k
                for l in j..(-1)..Array.max [|i+3; j-IL_MAX_SIZE; k+3|] do // right loop l, j
                    let tmp = model.internalLoop i k l j + V (k+1) (l-1)
                    if tmp > max then max <- tmp
            if j-i+1 > 3 then // try multiloop if there is enough room
                max <- Math.Max(max, IM i j)
            dpW.[i,j] <- max
            max
    and V i j = // i and j close a stem externally
        if i >= j then Double.MinValue
        elif dpV.[i,j] > Double.MinValue then dpV.[i,j]
        else
            let mutable max = Double.MinValue
            for k in 1..(j-i+1)/2 do // k is stem size
                let tmp = (model.stem i j k) + 
                            if i+k <= j-k then W (i+k) (j-k) else 0.0 // stems with no hp loop
                if tmp > max then max <- tmp
            dpV.[i,j] <- max
            max
    and IM i j = // there is an internal multiloop starting at i and ending at j
        if i > j then 0.0
        elif dpIM.[i,j] > System.Double.MinValue then dpIM.[i,j]
        else
            let mutable max = System.Double.MinValue
            for k in i..j do // k is start of dangle or stem
                let tmp =  IM i (k-1) + System.Math.Max(V k j + model.interalStemBonus k j, model.internalDangle k j)
                if tmp > max then max <- tmp
            dpIM.[i,j] <- max
            max
    and EM j = // external loop ending at j
        if j < 0 then 0.0
        elif dpEM.[j] > Double.MinValue then dpEM.[j]
        else
            let mutable max = Double.MinValue
            for k in 0..j do // k is start of dangle or stem
                let tmp = EM (k-1) + System.Math.Max(V k j + model.externalStemBonus k j, model.extenalDangle k j)
                if tmp > max then max <- tmp
            dpEM.[j] <- max
            max
    let tmp = EM (rna.Length-1)
    printfn "%A" dpEM
    printfn "%A" dpIM
    printfn "%A" dpV
    printfn "%A" dpW
    tmp