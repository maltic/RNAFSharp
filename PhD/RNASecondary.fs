module RNASecondary

/// A base pairing.
type Pair = 
    | GU = 0
    | UG = 1
    | GC = 2
    | CG = 3
    | AU = 4
    | UA = 5

/// All valid pairs.
let pairs = [| Pair.GU; Pair.UG; Pair.GC; Pair.CG; Pair.AU; Pair.UA |]

/// Check if a and b form a canonical base pair
let validPair a b = 
    let arr = [| a; b |] |> Array.sort
    match arr with
    | [| RNAPrimary.Base.A; RNAPrimary.Base.U |] -> true
    | [| RNAPrimary.Base.U; RNAPrimary.Base.G |] -> true
    | [| RNAPrimary.Base.G; RNAPrimary.Base.C |] -> true
    | _ -> false

/// Given two bases that constitute a valid base pair, returns the Pair representation.
let basesToPair a b = 
    match a, b with
    | RNAPrimary.Base.G, RNAPrimary.Base.U -> Pair.GU
    | RNAPrimary.Base.U, RNAPrimary.Base.G -> Pair.UG
    | RNAPrimary.Base.G, RNAPrimary.Base.C -> Pair.GC
    | RNAPrimary.Base.C, RNAPrimary.Base.G -> Pair.CG
    | RNAPrimary.Base.A, RNAPrimary.Base.U -> Pair.AU
    | RNAPrimary.Base.U, RNAPrimary.Base.A -> Pair.UA
    | _ -> failwith "Invalid pair"

/// A structural surface. Basically a chord in the circle graph representation of an RNA secondary structure.
type Surface = 
    | Reducible of Surface list * int * int * int // stem (i, j, size)
    | Irreducible of int * int // loop (i, j)

/// Parse an irreducible surface
let parseIS (str:string) i = 
    Irreducible(i, 
        (Seq.skipWhile (fun i -> i < str.Length && str.[i] = '.') (seq { i..str.Length }) |> Seq.head)-1)

/// Parse the exterior loop of an RNA secondary structure into surfaces
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


/// Free energy model
type Model = 
    { 
        stem : int -> int -> int -> double;
        internalLoop : int -> int -> int -> int -> double;
        lbulge : int -> int -> double;
        rbulge : int -> int -> double;
        hairpin : int -> int -> double;
        internalDangle : int -> int -> double;
        externalDangle : int -> int -> double;
        externalStemBonus : int -> int -> double;
        interalStemBonus : int -> int -> double;
    }


/// Apply an energy model to an rna secondary structure
let scoreExternalLoop model loop =
    let rec scoreInternalSurface = function
        | Reducible(children, i, j, sz) ->
            model.stem i j sz + 
                match children with
                | [Irreducible(i, j); (Reducible(_) as r); Irreducible(k, l)] -> 
                    model.internalLoop i j k l + scoreInternalSurface r
                | [Irreducible(i, j); (Reducible(_) as r)] -> 
                    model.lbulge i j + scoreInternalSurface r
                | [(Reducible(_) as r); Irreducible(i, j)] -> 
                    model.rbulge i j + scoreInternalSurface r
                | [Irreducible(i, j)] -> model.hairpin i j
                | [] -> model.hairpin (i+sz) (j-sz)
                | l -> List.fold (fun s t -> match t with 
                                                | Irreducible(i, j) -> model.internalDangle i j + s
                                                | Reducible(_, i, j, _) as r -> 
                                                    model.interalStemBonus i j + 
                                                        scoreInternalSurface r + s) 0.0 l
        | Irreducible(_) -> failwith "Impossible"
    List.fold (fun s t -> match t with 
                            | Irreducible(i, j) -> model.externalDangle i j + s
                            | Reducible(_, i, j, _) as r -> 
                                model.externalStemBonus i j + 
                                    scoreInternalSurface r + s) 0.0 loop

/// Enumerates all valid secondary structures                      
let validSecondaryStructures (rna:RNAPrimary.Base[]) = 
    let rec enumerate i stack (s:char[]) = 
        seq {
            if i = rna.Length then 
                if List.isEmpty stack then yield System.String(s)
            else
                match stack with
                | h :: t when validPair rna.[h] rna.[i] ->
                    s.[i] <- ')'
                    yield! enumerate (i+1) t s
                | _ -> ()
                s.[i] <- '('
                yield! enumerate (i+1) (i::stack) s
                s.[i] <- '.'
                yield! enumerate (i+1) stack s
        }
    enumerate 0 [] (Array.zeroCreate rna.Length)
