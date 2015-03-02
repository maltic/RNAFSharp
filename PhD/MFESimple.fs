module MFESimple

open System

type Parameterization = 
    { rna : RNAPrimary.Base []
      p : TurnerSimple.Parameters }
    
    member this.loopEnergy i k l j = 
        if i = k && l = j then Double.MaxValue
        elif k = l then this.p.hp.score this.rna i j
        elif i = k then this.p.bulge.score (j - l + 1)
        elif l = j then this.p.bulge.score (k - i + 1)
        else this.p.intern.score this.rna i k l j
    
    member this.W i j = 
        if i > j then Double.MaxValue
        else 
            Seq.min (seq { 
                         yield this.M i j 0
                         for k in i..j - 1 do
                             for l in k..j do
                                 yield this.loopEnergy i k l j + this.V k l
                     })
    
    member this.V i j = 
        if i >= j then Double.MaxValue
        else 
            Seq.min (seq { 
                         yield this.p.stack.scoreSurface this.rna i j + this.V (i + 1) (j - 1)
                         yield this.W (i + 1) (j - 1) + this.p.stack.init
                     })
    
    member this.M i j b = 
        if i > j then 
            if b = 2 then this.p.multi.b + this.p.multi.a
            else Double.MaxValue
        else 
            Seq.min (seq { 
                         yield this.p.multi.b + this.M i (j - 1) b
                         for k in i..j - 1 do
                             yield this.p.multi.c + this.M i (k - 1) (Math.Max(2, b + 1)) + this.V k j
                     })
