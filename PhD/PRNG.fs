module PRNG
open FSharpx.Collections

let stream seed generator = 
    let r = System.Random(seed)
    LazyList.ofSeq (Seq.initInfinite (fun _ -> generator r))

let doubles seed = stream seed (fun r -> r.NextDouble())

let ints seed = stream seed (fun r -> r.Next())