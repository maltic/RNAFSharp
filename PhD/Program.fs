
open System

[<EntryPoint>]
let main argv = 
    MFESimple.test 8 10000
    System.Console.ReadLine() |> ignore
    0 // return an integer exit code
