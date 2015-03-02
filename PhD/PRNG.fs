module PRNG
open FSharpx.Collections

let stream seed generator = 
    let r = System.Random(seed)
    LazyList.ofSeq (Seq.initInfinite (fun _ -> generator r))