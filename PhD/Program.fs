

type RNASurface = 
    | Reducible of RNASurface list * int * int * int // stem
    | Irreducible of int * int // loop

let parseIS (str:string) i = 
    Irreducible(i, 
        (Seq.skipWhile (fun i -> i < str.Length && str.[i] = '.') (seq { i..str.Length }) |> Seq.head)-1)

let parseSurfaces (str:string) = 
    let rec parse stack i =
        if i >= str.Length then // base case, must reach the end of the outer surface 
            match stack with
            | (s, elems) :: [] -> List.rev elems
            | _ -> failwith "Mismatched parens"
        else
            match str.[i], stack with
            | '.', (hi, h) :: t -> // irreducible surface
                match parseIS str i with
                | Irreducible(_, e) as s -> parse ((hi, s::h) :: t) (e+1)
                | _ -> failwith "Impossible"
            | '(', stack -> parse ((i, []) :: stack) (i+1) // push reducible surface on stack
            | ')', (hi, h) :: (mi, m) :: t -> // pop reducible from stack, merge with next stack element (the parent)
                match h with
                | [Reducible(children, _, _, c)] -> 
                    parse ((mi, (Reducible(children, hi, i, c+1) :: m)) :: t) (i+1) // continues a stem
                | _ -> 
                    parse ((mi, (Reducible(List.rev h, hi, i, 1) :: m)) :: t) (i+1) // closes a loop
            | ')' , _ :: [] -> failwith "Mismatched parens" // should never need to pop the outer surface
            | _ -> failwith "Invalid character"
    parse [(-1, [])] 0

let stemScore i j k = 1.0
let internalLoopScore i k j l = 1.0
let bulgeScore i j = 1.0
let hairpinLoopScore i j = 1.0
let internalDangleScore i j = 1.0
let extenalDangleScore i kj= 1.0
let externalStemBonus i j = 1.0
let interalStemBonus i j = 1.0

type Model = { 
    stem : int -> int -> int -> double;
    internalLoop : int -> int -> int -> int -> double;
    lbulge : int -> int -> double;
    rbulge : int -> int -> double;
    hairpin : int -> int -> double;
    internalDangle : int -> int -> double;
    extenalDangle : int -> int -> double;
    externalStemBonus : int -> int -> double;
    interalStemBonus : int -> int -> double;
}

let scoreExternalLoop model loop =
    let rec scoreInternalSurface model = function
        | Reducible(children, i, j, sz) ->
            model.stem i j sz + 
                match children with
                | [Irreducible(i, j); (Reducible(_) as r); Irreducible(k, l)] -> 
                    model.internalLoop i j k l + scoreInternalSurface model r
                | [Irreducible(i, j); (Reducible(_) as r)] -> 
                    model.lbulge i j + scoreInternalSurface model r
                | [(Reducible(_) as r); Irreducible(i, j)] -> 
                    model.rbulge i j + scoreInternalSurface model r
                | [Irreducible(i, j)] -> model.hairpin i j
                | l -> List.fold (fun s t -> match t with 
                                                | Irreducible(i, j) -> model.internalDangle i j + s
                                                | r -> scoreInternalSurface model r + s) 0.0 l
        | Irreducible(_) -> failwith "Impossible"
    List.fold (fun s t -> match t with 
                            | Irreducible(i, j) -> model.extenalDangle i (j-i+1) + s
                            | r -> scoreInternalSurface model r + s) 0.0 loop

let zuker (rna:string) model = 
    let dpW = Array2D.create rna.Length rna.Length System.Double.MinValue
    let dpV = Array2D.create rna.Length rna.Length System.Double.MinValue
    let dpIM = Array2D.create rna.Length rna.Length System.Double.MinValue
    let dpEM = Array.create rna.Length System.Double.MinValue
    let IL_MAX_SIZE = 30
    let rec W i j = // i and j begin a loop (i-1 and j+1 are bonded)
        if i > j then 0.0
        elif dpW.[i,j] > System.Double.MinValue then dpW.[i,j]
        else
            let mutable max = model.hairpin i j // hairpin
            for k in i..(j-1) do // k is the start of the succeeding stem, bulge from i, k-1
                let tmp = model.lbulge i (k-1) + V k j
                if tmp > max then max <- tmp
            for k in j..(-1)..i+1 do // k is the end of the precding stem, bulge from k+1, j
                let tmp = model.rbulge (k+1) j + V i k
                if tmp > max then max <- tmp
            // internal loops
            for k in i..System.Math.Min(j-3, i+IL_MAX_SIZE) do // left loop i, k
                for l in j..(-1)..Array.max [|i+3; j-IL_MAX_SIZE; k+3|] do // right loop l, j
                    let tmp = model.internalLoop i k l j + V (k+1) (l-1)
                    if tmp > max then max <- tmp
            dpW.[i,j] <- System.Math.Max(max, IM i j) // try multiloop
            dpW.[i,j]
    and V i j = // i and j close a stem externally
        if i >= j then 0.0
        elif dpV.[i,j] > System.Double.MinValue then dpV.[i,j]
        else
            let mutable max = System.Double.MinValue
            for k in 1..(j-i+1)/2 do // k is stem size
                let tmp = (model.stem i j k) + W (i+k) (j-k)
                if tmp > max then max <- tmp
            dpV.[i,j] <- max
            max
    and IM i j = // there is an internal multiloop starting at i and ending at j
        if i > j then 0.0
        elif dpIM.[i,j] > System.Double.MinValue then dpIM.[i,j]
        else
            let mutable max = System.Double.MinValue
            for k in i..j do // k is start of dangle or stem
                let tmp =  IM i (k-1) + System.Math.Max(V k j + model.interalStemBonus k j, model.internalDangle k j)
                if tmp > max then max <- tmp
            dpIM.[i,j] <- max
            max
    and EM j = // external loop ending at j
        if j < 0 then 0.0
        elif dpEM.[j] > System.Double.MinValue then dpEM.[j]
        else
            let mutable max = System.Double.MinValue
            for k in 0..j do // k is start of dangle or stem
                let tmp = EM (k-1) + System.Math.Max(V k j + model.externalStemBonus k j, model.extenalDangle k j)
                if tmp > max then max <- tmp
            dpEM.[j] <- max
            max
    EM (rna.Length-1)

type GeneticAlgorithm<'a> = { 
    breeder : System.Random -> 'a -> 'a -> 'a; 
    mutator: System.Random -> 'a -> 'a;
    creator : System.Random -> 'a;
    fitness : 'a -> float;
    elitism : float;
    newBlood : float;
    spawn : float;
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

let testBreeder _ a b = (a+b)/2.0
let testMutate (r:System.Random) i = i + (r.NextDouble() - 0.5)
let testFitness i = System.Math.Abs(12.32 - i)
let testCreate (r:System.Random) = r.NextDouble()
let testGA = {breeder = testBreeder; mutator = testMutate; fitness = testFitness;
    elitism = 0.1; mutate = 0.8; spawn = 0.0; selection = 0.4; newBlood = 0.05;
    creator = testCreate}
     
let evolve ga popSize steps = 
    let r = System.Random()
    Seq.fold (fun (s:'a[]) _ -> printfn "%f" (ga.fitness s.[0]); runGA r ga s) 
        (Array.init popSize (fun _ -> ga.creator r)) 
        (seq {1..steps})

[<EntryPoint>]
let main argv = 
    let m = {
        stem=stemScore; internalDangle = internalDangleScore; extenalDangle = extenalDangleScore; lbulge = bulgeScore;
        internalLoop = internalLoopScore; hairpin=hairpinLoopScore; rbulge = bulgeScore; externalStemBonus = externalStemBonus; interalStemBonus = interalStemBonus} 
    //System.Console.ReadLine() |> parseSurfaces |> scoreExternalLoop m |> printfn "%A"
    //evolve testGA 100 300 |> printfn "%A"
    //let n = 10000000
    //printfn "%A" (parseSurfaces (new string [| for i in 1..n do yield if i <= (n/2) then '(' else ')'|]))
    let sw = System.Diagnostics.Stopwatch()
    sw.Start()
    zuker (new string (Array.create 251 ' ')) m |> printfn "%A"
    sw.Stop()
    printfn "%d" sw.ElapsedMilliseconds
    let sw = System.Diagnostics.Stopwatch()
    sw.Start()
    zuker (new string (Array.create 251 ' ')) m |> printfn "%A"
    sw.Stop()
    printfn "%d" sw.ElapsedMilliseconds
    let sw = System.Diagnostics.Stopwatch()
    sw.Start()
    zuker (new string (Array.create 251 ' ')) m |> printfn "%A"
    sw.Stop()
    printfn "%d" sw.ElapsedMilliseconds
    System.Console.Read() |> ignore
    0 // return an integer exit code
