

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
let bulgeScore i k = 1.0
let hairpinLoopScore i k = 1.0
let internalDangleScore i k = 1.0
let extenalDangleScore i k = 1.0

type Model = { 
    stem : int -> int -> int -> double;
    internalLoop : int -> int -> int -> int -> double;
    lbulge : int -> int -> double;
    rbulge : int -> int -> double;
    hairpin : int -> int -> double;
    internalDangle : int -> int -> double;
    extenalDangle : int -> int -> double;
}

let scoreExternalLoop model loop =
    let rec scoreInternalSurface model = function
        | Reducible(children, i, j, sz) ->
            model.stem i j sz + 
                match children with
                | [Irreducible(i, j); (Reducible(_) as r); Irreducible(k, l)] -> 
                    model.internalLoop i (j-i+1) k (l-k+1) + scoreInternalSurface model r
                | [Irreducible(i, j); (Reducible(_) as r)] -> 
                    model.lbulge i (j-i+1) + scoreInternalSurface model r
                | [(Reducible(_) as r); Irreducible(i, j)] -> 
                    model.rbulge i (j-i+1) + scoreInternalSurface model r
                | [Irreducible(i, j)] -> model.hairpin i (j-i+1)
                | l -> List.fold (fun s t -> match t with 
                                                | Irreducible(i, j) -> model.internalDangle i (j-i+1) + s
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
    let rec W i j = // i and j begin a loop externally (i-1 and j+1 are bonded)
        if i > j then 0.0
        elif dpW.[i,j] > System.Double.MinValue then dpW.[i,j]
        else
            let sz = j-i+1
            let mutable max = model.hairpin i sz
            // bulges
            for k in 1..sz do
                let tmp = System.Math.Max(model.lbulge i k + V (i+k) j, 
                            model.rbulge (j-k+1) k + V i (j-k))
                if tmp > max then max <- tmp
            // internal loops
            let maxILSZ = System.Math.Max(sz, 30)
            for k in 1..maxILSZ do
                for l in 1..maxILSZ do
                    if i+k-1 < j-l+1 then
                        let tmp = model.internalLoop i k j l + V (i+k) (j-l)
                        if tmp > max then max <- tmp
            // try multiloop
            dpW.[i,j] <- System.Math.Max(max, IM i j)
            dpW.[i,j]
    and V i j = // i and j close a stem externally
        if i >= j then 0.0
        elif dpV.[i,j] > System.Double.MinValue then dpV.[i,j]
        else
            let mutable max = System.Double.MinValue
            for k in 1..(j-i-1) do
                if i+k-1 < j-k+1 then
                    let tmp = (model.stem i j k) + W (i+k) (j-k)
                    if tmp > max then max <- tmp
            dpV.[i,j] <- max
            max
    and IM i j = // there is an internal multiloop starting at i and ending at j
        if i > j then 0.0
        elif dpIM.[i,j] > System.Double.MinValue then dpIM.[i,j]
        else
            let mutable max = System.Double.MinValue
            for k in 1..(j-i+1) do
                let tmp =  System.Math.Max(V (j-k+1) j + IM i (j-k), 
                            model.internalDangle (j-k+1) k + IM i (j-k))
                if tmp > max then max <- tmp
            dpIM.[i,j] <- max
            max
    and EM j = // external loop ending at j
        if j < 0 then 0.0
        elif dpEM.[j] > System.Double.MinValue then dpEM.[j]
        else
            let mutable max = System.Double.MinValue
            for k in 0..j do
                let tmp = EM (k-1) + System.Math.Max(V k j, model.extenalDangle k (j-k+1))
                if tmp > max then max <- tmp
            dpEM.[j] <- max
            max
    EM (rna.Length-1)

[<EntryPoint>]
let main argv = 
    let m = {
        stem=stemScore; internalDangle = internalDangleScore; extenalDangle = extenalDangleScore; lbulge = bulgeScore;
        internalLoop = internalLoopScore; hairpin=hairpinLoopScore; rbulge = bulgeScore;} 
    System.Console.ReadLine() |> parseSurfaces |> scoreExternalLoop m |> printfn "%A"
    System.Console.Read() |> ignore
    //let n = 10000000
    //printfn "%A" (parseSurfaces (new string [| for i in 1..n do yield if i <= (n/2) then '(' else ')'|]))
    0 // return an integer exit code
