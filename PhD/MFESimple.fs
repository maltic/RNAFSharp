module MFESimple

open System

type FETables(rna : RNAPrimary.Base [], p : TurnerSimple.Parameters)  =
    
    let sz = rna.Length
    let Vt = Array2D.create sz sz Double.NaN
    let Mt = Array3D.create sz sz 3 Double.NaN
    let Et = Array.create sz Double.NaN
    
    let loopEnergy i k l j = 
        if i + 1 = j then 0.0 // empty hairpin loops
        elif k = l then p.hp.score rna i j
        elif not (RNASecondary.validPair (rna.[k]) (rna.[l])) then Double.MaxValue
        elif k - 1 = i && l + 1 = j then p.stack.scoreSurface rna i j
        elif k - 1 = i then p.bulge.score (j - l - 1) + p.stack.init
        elif l + 1 = j then p.bulge.score (k - i - 1) + p.stack.init
        else p.intern.score rna i k l j + p.stack.init
    
    let rec V i j = 
        if i >= j then Double.MaxValue
        elif not (RNASecondary.validPair (rna.[i]) (rna.[j])) then Double.MaxValue
        elif not (Double.IsNaN Vt.[i,j]) then Vt.[i,j]
        else 
            Vt.[i,j] <- Seq.min (seq { 
                         for k in i + 1..j - 2 do
                             for l in k + 1..j-1 do
                                 yield loopEnergy i k l j + V k l
                         yield loopEnergy i -1 -1 j
                         yield M (i + 1) (j - 1) 0
                     })
            Vt.[i,j]
    
    and M i j b = 
        if i > j then 
            if b = 2 then p.multi.a
            else Double.MaxValue
        elif not (Double.IsNaN Mt.[i,j,b]) then Mt.[i,j,b]
        else 
            Mt.[i,j,b] <- Seq.min (seq { 
                         yield p.multi.b + M i (j - 1) b
                         for k in i..j - 1 do
                             yield M i (k - 1) (Math.Min(2, b + 1)) + V k j + p.multi.c 
                                   + p.stack.init
                     })
            Mt.[i,j,b]
    
    and E i = 
        if i < 0 then p.multi.a
        elif not (Double.IsNaN Et.[i]) then Et.[i]
        else 
            Et.[i] <- Seq.min (seq { 
                            yield p.multi.b + E (i - 1)
                            for j in 0..i - 1 do
                                yield E (j - 1) + p.multi.c + V j i + p.stack.init
                        } (*|> fun s -> printfn "%A %d" s i; s*))
            Et.[i]

    do E (rna.Length-1) |> ignore

    member this.ETable = Et
    member this.VTable = Vt
    member this.MTable = Mt

let test sz iterations = 
    let r = System.Random()
    
    let randomSequence sz = 
        let b = "augc"
        new string [| for i in 1..sz -> b.[r.Next(3)] |]
    
    let testSequence s = 
        let rna = 
            s
            |> RNAPrimary.parse
            |> Seq.toArray
        let p = PRNG.stream DateTime.Now.Millisecond (fun r -> r.NextDouble() * 10.0 - 5.0) |> TurnerSimple.Parameters.ofSeq
        let par = FETables(rna, p)
        
        let ssScore = RNASecondary.parseSurfaces >> p.score rna
        
        let best = 
            rna
            |> RNASecondary.validSecondaryStructures
            |> Seq.minBy ssScore
        
        let bests = ssScore best
        let zuks = par.ETable.[(rna.Length-1)]
        if System.Math.Abs(bests - zuks) > 0.00001 then 
            printfn "%A vs %A \n\n %A \n\n %A \n\n %A" bests zuks best rna p
            printfn "%A" par.ETable
            printfn "%A" par.VTable
            printfn "%A" par.MTable
            printfn "%A" (FETables(rna, p))
    
    for i in 1..iterations do
        testSequence (randomSequence sz)
