module MFE

open RNASecondary
open System

let zuker (rna:RNAPrimary.Base[]) model = 
    let dpW = Array2D.create rna.Length rna.Length Double.MinValue
    let dpV = Array2D.create rna.Length rna.Length Double.MinValue
    let dpIM = Array4D.create rna.Length rna.Length 3 2 Double.MinValue
    let dpEM = Array2D.create rna.Length 2 Double.MinValue
    let IL_MAX_SIZE = 30
    let rec W i j = // i and j begin a loop (i-1 and j+1 are bonded)
        if dpW.[i,j] > Double.MinValue then dpW.[i,j]
        else
            let mutable max = model.hairpin i j // hairpin, lack of base-case allows stems without hp loops
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
                max <- Math.Max(max,    // try starting with/without stem
                    Math.Max(IM i j 2 0 |> float, IM i j 2 1))
            dpW.[i,j] <- max
            max
    and V i j = // i and j close a stem externally
        if i >= j then Double.MinValue
        elif dpV.[i,j] > Double.MinValue then dpV.[i,j]
        else
            let mutable max = Double.MinValue
            for k in 1..(j-i+1)/2 do // k is stem size
                let tmp = (model.stem i j k) + W (i+k) (j-k)
                if tmp > max then max <- tmp
            dpV.[i,j] <- max
            max
    /// there is an internal multiloop starting at i and ending at j, rem is the remaining stems to place
    and IM i j rem stem =
        if i > j then if rem = 0 then 0.0 else Double.MinValue
        elif dpIM.[i,j,rem,stem] > System.Double.MinValue then dpIM.[i,j,rem,stem]
        else
            let mutable max = System.Double.MinValue
            for k in i..j do // k is start of dangle or stem
                let tmp = 
                    if stem = 1 then
                        let nextRem = Math.Max(0, rem-1)
                        V k j + model.interalStemBonus k j //coaxial and non coaxial
                            + Math.Max(IM i (k-1) nextRem 0, IM i (k-1) nextRem 1)
                    else
                        IM i (k-1) rem 1 + model.internalDangle k j
                if tmp > max then max <- tmp
            dpIM.[i,j,rem,stem] <- max
            max
    and EM j stem = // external loop ending at j
        if j < 0 then 0.0
        elif dpEM.[j,stem] > Double.MinValue then dpEM.[j,stem]
        else
            let mutable max = Double.MinValue
            for k in 0..j do // k is start of dangle or stem
                let tmp =
                    if stem = 1 then 
                        V k j + model.externalStemBonus k j + 
                            Math.Max(EM (k-1) 0, EM (k-1) 1) // try coaxial stack, and stack-stem
                    else 
                        model.externalDangle k j + EM (k-1) 1
                if tmp > max then max <- tmp
            dpEM.[j,stem] <- max
            max
    Math.Max(EM (rna.Length-1) 0, EM (rna.Length-1) 1)