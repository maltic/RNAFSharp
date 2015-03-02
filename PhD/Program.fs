

[<EntryPoint>]
let main argv = 
    MFESimple.test 4 100
    System.Console.ReadLine() |> ignore
    0 // return an integer exit code
