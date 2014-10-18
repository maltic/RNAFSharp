
type RNASurface = 
    | Reducible of RNASurface list * int * int
    | Irreducible of int * int

let parseIS (str:string) i = 
    Irreducible(i, 
        (Seq.skipWhile (fun i -> i < str.Length && str.[i] = '.') (seq { i..str.Length }) |> Seq.head)-1)

let parseSurfaces (str:string) = 
    let rec parse stack i =
        if i >= str.Length then // base case, must reach the end of the outer surface 
            match stack with
            | (s, elems) :: [] -> Reducible(List.rev elems, s, i)
            | _ -> failwith "Mismatched parens"
        else
            match str.[i], stack with
            | '.', (hi, h) :: t -> // irreducible surface
                match parseIS str i with
                | Irreducible(_, e) as s -> parse ((hi, s::h) :: t) (e+1)
                | _ -> failwith "Impossible"
            | '(', stack -> parse ((i, []) :: stack) (i+1) // push reducible surface on stack
            | ')', (hi, h) :: (mi, m) :: t -> // pop reducible from stack, merge with next stack element (the parent)
                parse ((mi, (Reducible(List.rev h, hi, i) :: m)) :: t) (i+1)
            | ')' , _ :: [] -> failwith "Mismatched parens" // should never need to pop the outer surface
            | _ -> failwith "Invalid character"
    parse [(-1, [])] 0

let stemScore i j k = k
let internalLoopScore i k j l = k + l
let bulgeScore i k = k
let hairpinLoopScore i k = k
let extenalDangleScore i k = k

let zuker (rna:string) = 
    let dpW = Array2D.create rna.Length rna.Length -1
    let dpV = Array2D.create rna.Length rna.Length -1
    let dpIM = Array2D.create rna.Length rna.Length -1
    let rec W i j = 
        if i > j then 0
        elif dpW.[i,j] > -1 then dpW.[i,j]
        else
            let sz = j-i+1
            let mutable max = hairpinLoopScore i sz
            // bulges
            for k in 1..sz do
                let tmp = System.Math.Max(bulgeScore i k + V (i+k) j, 
                            bulgeScore (j-k+1) k + V i (j-k))
                if tmp > max then max <- tmp
            let maxILSZ = System.Math.Max(sz, 30)
            for k in 1..maxILSZ do
                for l in 1..maxILSZ do
                    if i+k-1 < j-l+1 then
                        let tmp = internalLoopScore i k j l + V (i+k) (j-l)
                        if tmp > max then max <- tmp
            // try multiloop
            dpW.[i,j] <- System.Math.Max(max, IM i j)
            dpW.[i,j]
    and V i j = 
        if i >= j then 0
        elif dpV.[i,j] > -1 then dpV.[i,j]
        else
            let mutable max = -1
            for k in 1..(j-i) do
                if i+k-1 < j-k+1 then
                    let tmp = (stemScore i j k) + W (i+k) (j-k)
                    if tmp > max then max <- tmp
            dpV.[i,j] <- max
            max
    and IM i j = 
        if i > j then 0
        elif dpIM.[i,j] > -1 then dpIM.[i,j]
        else
            let mutable max = -1
            for k in 1..(j-i+1) do
                let tmp =  System.Math.Max(V (j-k+1) j + IM i (j-k), 
                            extenalDangleScore (j-k+1) k + IM i (j-k))
                if tmp > max then max <- tmp
            dpIM.[i,j] <- max
            max
    IM 0 (rna.Length-1)

[<EntryPoint>]
let main argv = 
    let sw = System.Diagnostics.Stopwatch()
    sw.Start()
    zuker (new string [| for i in 1..110 do yield ' ' |]) |> printfn "%A"
    sw.Stop()
    printfn "%A" (sw.ElapsedMilliseconds)
    System.Console.Read() |> ignore
    //let n = 10000000
    //printfn "%A" (parseSurfaces (new string [| for i in 1..n do yield if i <= (n/2) then '(' else ')'|]))
    0 // return an integer exit code
