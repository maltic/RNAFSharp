module RNASecondary


type Surface = 
    | Reducible of Surface list * int * int * int // stem
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

