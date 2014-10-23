module RNAPrimary

let BASES = 4

type Base = A = 0 | U = 1 | G = 2 | C = 3

/// Check if a and b form a canonical base pair
let validPair a b = 
    let arr = [| a; b |] |> Array.sort
    match arr with
    | [| Base.A; Base.U |] -> true
    | [| Base.U; Base.G |] -> true
    | [| Base.G; Base.C |] -> true
    | _ -> false

/// Parse a string into an rna sequence, reports errors
let parse (str:string) =
    seq { for c in str ->
                match System.Char.ToUpper(c) with
                | 'A' -> Base.A
                | 'U' -> Base.U
                | 'G' -> Base.G
                | 'C' -> Base.C
                | _ -> failwith "Invalid nucleotide" }