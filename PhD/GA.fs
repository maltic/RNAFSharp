module GA
open FSharpx.Collections
open System

type Configuration<'a> = { 
    breeder : LazyList<float> -> 'a -> 'a -> 'a; 
    mutator: LazyList<float> -> 'a -> 'a;
    creator : LazyList<float> -> 'a;
    fitness : 'a -> float;
    elitism : float;
    newBlood : float;
    mutate : float;
    selection : float;
}

let run randStream seed config population = 
    let popSize = Array.length population
    let elitism = Math.Round(float popSize * config.elitism) |> int
    let selected = Math.Round (config.selection * float popSize) |> int
    let newBlood = Math.Round (config.newBlood * float popSize) |> int
    seq {
        for generation in Seq.initInfinite (fun i -> i) do
            let next = 
                Array.Parallel.init popSize 
                    (fun i ->
                        let r = PRNG.stream (seed+generation+i) (fun (r:Random) -> r.NextDouble())
                        if i < elitism then population.[i]
                        elif i < elitism + newBlood then config.creator r
                        else
                            if LazyList.head r > config.mutate then
                                config.breeder (LazyList.tail r) population.[i % selected] population.[i % selected]
                            else
                                config.mutator (LazyList.tail r) population.[i % selected]
                    )
            Array.sortInPlaceBy config.fitness next
            yield next
    }

// some helper functions
let maybeMutate rate range r a = 
        if LazyList.head r > rate then
            a
        else
            a + ((Seq.nth 1 r * range) - range/2.0)

let maybeTake rate r a b = 
    if LazyList.head r > rate then
        a
    else
        b