﻿

[<EntryPoint>]
let main argv = 
    MFESimple.test 1 2
    System.Console.ReadLine() |> ignore
    0 // return an integer exit code
