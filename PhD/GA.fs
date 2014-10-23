module GA


type GeneticAlgorithm<'a> = { 
    breeder : System.Random -> 'a -> 'a -> 'a; 
    mutator: System.Random -> 'a -> 'a;
    creator : System.Random -> 'a;
    fitness : 'a -> float;
    elitism : float;
    newBlood : float;
    mutate : float;
    selection : float;
}

let runGA r ga population = 
    let popSize = Array.length population
    let elitism = System.Math.Round(float popSize * ga.elitism) |> int
    let selected = System.Math.Round (ga.selection * float popSize) |> int
    let newBlood = System.Math.Round (ga.newBlood * float popSize) |> int
    let next = 
        Array.init popSize 
            (fun i ->
                if i < elitism then population.[i]
                elif i < elitism + newBlood then ga.creator r
                else
                    if r.NextDouble() > ga.mutate then
                        ga.breeder r population.[i % selected] population.[i % selected]
                    else
                        ga.mutator r population.[i % selected]
            )
    Array.sortInPlaceBy ga.fitness next
    next