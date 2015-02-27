module RNAPrimary

type Base = A = 0 | U = 1 | G = 2 | C = 3

// TODO This should be made immutable
let bases = [| Base.A; Base.U; Base.G; Base.C |]

/// Parse a string into an rna sequence, reports errors
let parse (str:string) =
    seq { for c in str ->
                match System.Char.ToUpper(c) with
                | 'A' -> Base.A
                | 'U' -> Base.U
                | 'G' -> Base.G
                | 'C' -> Base.C
                | _ -> failwith "Invalid nucleotide" }