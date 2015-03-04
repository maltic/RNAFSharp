/// Tools for doing priority first free energy searching
module PFFE
open FSharpx.Collections

let search feScorer rna =    
    let stems = RNASecondary.Stem.enumerate rna |> Array.ofSeq
    let heap = Heap.ofSeq false (Array.map (fun stem -> (0, stem)) stems)
    Seq.unfold (fun pq -> if Heap.isEmpty pq then None else Some((), pq)) heap

