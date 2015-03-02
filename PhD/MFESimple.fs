module MFESimple

let rec W e i j =
    0.0
and V rna (e:TurnerSimple.Parameters) i j =
    Seq.min (seq {
        yield e.stack.score rna i j + V rna e (i+1) (j-1)
        yield W e i j
    })
and M e i j = 
    0.0