module MFESimple

open System

type Parameterization = 
    { rna : RNAPrimary.Base []
      p : TurnerSimple.Parameters }
    
    member this.loopEnergy i k l j = 
        if not (RNASecondary.validPair (this.rna.[i]) (this.rna.[j])) then Double.MaxValue
        elif i + 1 = j then 0.0 // empty hairpin loops
        elif k = l then this.p.hp.score this.rna i j
        elif not (RNASecondary.validPair (this.rna.[k]) (this.rna.[l])) then Double.MaxValue
        elif k - 1 = i && l + 1 = j then this.p.stack.scoreSurface this.rna i j
        elif k = i + 1 then this.p.bulge.score (j - l - 1)
        elif l = j - 1 then this.p.bulge.score (k - i - 1)
        else this.p.intern.score this.rna i k l j
    
    member this.V i j = 
        if i >= j then Double.MaxValue
        else 
            Seq.min (seq { 
                         for k in i + 1..j - 2 do
                             for l in k + 1..j - 1 do
                                 yield this.loopEnergy i k l j + this.V k l
                         yield this.loopEnergy i i i j
                         yield this.M (i + 1) (j - 1) 0
                     })
    
    member this.M i j b = 
        if i > j then 
            if b = 2 then this.p.multi.a
            else Double.MaxValue
        else 
            Seq.min (seq { 
                         yield this.p.multi.b + this.M i (j - 1) b
                         for k in i..j - 1 do
                             yield this.M i (k - 1) (Math.Max(2, b + 1)) + this.V k j + this.p.multi.c 
                                   + this.p.stack.init
                     })
    
    member this.E i = 
        if i < 0 then this.p.multi.a
        else 
            Seq.min (seq { 
                         yield this.p.multi.b + this.E(i - 1)
                         for j in 0..i - 1 do
                             let a = this.V j i
                             let g = ()
                             yield this.E(j - 1) + this.p.multi.c + this.V j i + this.p.stack.init
                     })

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
        
        let par = 
            { rna = rna
              p = 
                  TurnerSimple.Parameters.ofSeq 
                      (PRNG.stream (System.DateTime.Now.Millisecond) (fun r -> r.NextDouble())) }
        
        let ssScore = RNASecondary.parseSurfaces >> par.p.score par.rna
        
        let best = 
            rna
            |> RNASecondary.validSecondaryStructures
            |> Seq.minBy ssScore
        
        let bests = ssScore best
        let zuks = par.E(rna.Length - 1)
        if System.Math.Abs(bests - zuks) > 0.00001 then 
            printfn "%A vs %A \n\n %A \n\n %A \n\n %A" bests zuks best rna par
    
    for i in 1..iterations do
        testSequence (randomSequence sz)
