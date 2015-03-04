
open System

[<EntryPoint>]
let main argv = 
    PFFE.test 40
    //MFESimple.test 10 1000
    //System.Console.ReadLine() |> RNAPrimary.parse |> Array.ofSeq |> RNASecondary.Stem.enumerate |> List.ofSeq |> printfn "%A"
    //System.Console.ReadLine() |> ignore
    0 // return an integer exit code
