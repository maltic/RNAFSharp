
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



[<EntryPoint>]
let main argv = 
    let n = 10000000
    printfn "%A" (parseSurfaces (new string [| for i in 1..n do yield if i <= (n/2) then '(' else ')'|]))
    0 // return an integer exit code
