/// Tools for doing priority first free energy searching
module PFFE

open FSharpx.Collections
open RNASecondary

let search feScorer rna = 
    let stems = Stem.enumerate rna |> Array.ofSeq
    let score = RNASecondary.parseSurfaces >> feScorer
    let stemSS stem = 
        let i, j, k = stem.i, stem.j, stem.pairs
        new string [| for n in 0..rna.Length-1 -> if n >= i && n < i+k then '(' elif n > j-k && n <= j then ')' else  '.' |]
    let mergeSSs (ssa:string) (ssb:string) = 
        if Seq.zip ssa ssb |> Seq.exists (fun (a,b) -> a <> '.' && b <> '.') then None
        // this needs to change so that it can enumerate multiloops
        elif (ssa.LastIndexOf('(') < ssb.IndexOf('(') && ssa.IndexOf(')') > ssb.LastIndexOf(')')) || (ssa.IndexOf('(') > ssb.LastIndexOf(')') || ssa.LastIndexOf(')') < ssb.IndexOf('(')) then
            Some (new string [| for i in 0..ssa.Length-1 -> if ssa.[i] <> '.' then ssa.[i] elif ssb.[i] <> '.' then ssb.[i] else '.'|])
        else
            None
    let stemSSs = Array.map stemSS stems
    let heap = 
        Heap.ofSeq false (Array.map (fun ss -> (score ss, ss)) stemSSs)
    Seq.unfold (fun (best, pq) -> 
        if Heap.isEmpty pq then None
        else
            let hdscore, hdstr = Heap.head pq
            let nextBest = if best < hdscore then best else hdscore
            let next = 
                stemSSs 
                |> Seq.map (mergeSSs hdstr) 
                |> Seq.filter Option.isSome 
                |> Seq.map (fun s -> (score (Option.get s), Option.get s)) 
                |> Heap.ofSeq false
            Some((nextBest, Heap.head pq), (nextBest, Heap.merge (Heap.tail pq) next) )) (System.Double.MaxValue, heap)

let test sz = 
    let r = System.Random()
    
    let randomSequence sz = 
        let b = "augc"
        new string [| for i in 1..sz -> b.[r.Next(3)] |]

    let rna = randomSequence sz |> RNAPrimary.parse |> Array.ofSeq
    let p = PRNG.stream System.DateTime.Now.Millisecond (fun r -> r.NextDouble() * 10.0 - 5.0) |> TurnerSimple.Parameters.ofSeq
    let fet = MFESimple.FETables(rna, p)
    printfn "Target: %f" (fet.ETable.[rna.Length-1])
    search (p.score rna) rna |> Seq.iter (printfn "%A")

